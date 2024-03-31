from __future__ import annotations

import argparse
import csv
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import Future as ComputeFuture
from typing import Sequence

from time import perf_counter

from proxystore.connectors.file import FileConnector
from proxystore.store import Store
from proxystore.store.future import ProxyFuture

from globus_compute_sdk import Executor

from genomes.individuals import processing_2
from genomes.individuals import processing_chrom_parts
from genomes.sifting import sifting
from genomes.mutation_overlap import run_moverlap
from genomes.frequency import run_frequency
from genomes.utils import Bench


def process_benchmarks(benchmarks):
    df = pd.DataFrame(benchmarks)
    colours = {
        "processing_chrom_parts": "b",
        "processing_2": "r",
        "sifting": "m",
        "mutation_overlap": "y",
        "frequency": "g",
    }
    df["colour"] = df["task"].map(colours)
    df["norm_start"] = df["start"] - df["start"].min()
    return df


class Workflow:
    populations_dir: str = "data/populations/"
    datafile: str = "data.csv"
    dataset: str = "20130502"
    columns: str = "columns.txt"
    ind_jobs: int
    c_nums: list
    individuals_files: list
    sifted_files: list
    output_fns: dict
    populations: list[str]
    debug: bool
    fraction: float

    def __init__(
        self,
        data_dir: str,
        results_dir: str,
        ind_jobs: int = 1,
        debug: bool = False,
        fraction: float = 1.0,
    ) -> None:
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.debug = debug
        self.fraction = fraction
        self.ind_jobs = ind_jobs

        populations_dir = os.path.join(data_dir, self.populations_dir)
        self.populations = [os.path.basename(f) for f in os.listdir(populations_dir)]
        self.c_nums = []
        self.individuals_files = []
        self.sifted_files = []
        self.output_fns = {}


def run_gc_workflow(cfg: Workflow, executor: Executor) -> list[Bench]:
    benchmarks: list[Bench] = []

    with open(cfg.datafile, "r") as f:
        for row in csv.reader(f):
            base_file = row[0]
            threshold = int(cfg.fraction * int(row[1]))
            # To ensure we do not create too many individuals jobs
            cfg.ind_jobs = min(cfg.ind_jobs, threshold)
            step = threshold // cfg.ind_jobs
            rest = threshold % cfg.ind_jobs
            if rest != 0:
                raise ValueError(
                    f"For file {base_file}: required individuals jobs "
                    f"{cfg.ind_jobs} does not divide the number of rows "
                    f"{threshold}."
                )

            counter = 1

            output_files = []

            # Individuals Jobs
            f_individuals = base_file
            base_file_path = os.path.join(
                cfg.data_dir, "data", cfg.dataset, f_individuals
            )

            c_num = base_file[base_file.find("chr") + 3 :]
            c_num = c_num[0 : c_num.find(".")]
            cfg.c_nums.append(c_num)

            if cfg.debug:
                cfg.c_nums = [1, 2]
                threshold = 2

            individuals_filename = "chr%sn" % c_num

            while counter < threshold:
                stop = counter + step

                out_name = os.path.join("chr%sn-%s-%s.tar.gz" % (c_num, counter, stop))
                output_files.append(out_name)

                f_chrom_parts = executor.submit(
                    processing_chrom_parts,
                    inputfile=base_file_path,
                    columnfile=os.path.join(
                        cfg.data_dir, "data", cfg.dataset, cfg.columns
                    ),
                    c=c_num,
                    counter=counter,
                    stop=stop,
                    total=threshold,
                    results_dir=cfg.results_dir,
                )

                chrom_parts = f_chrom_parts.result()
                chrom_bench = chrom_parts[0]
                benchmarks.append(chrom_bench)

                # print(f_chrom_parts)
                print("Done obtaining chromosome parts")
                future_chrn_df = []
                f_chrn_bench = []
                chrp_elements = chrom_parts[1][2]

                if cfg.debug:
                    chrp_elements = chrp_elements[0:5]

                futures = []
                for chrp in chrp_elements:
                    f_chrn_df = executor.submit(
                        processing_2,
                        chrp_element=chrp,
                        data=chrom_parts[1][0],
                        ndir=chrom_parts[1][1],
                        c=c_num,
                        results_dir=cfg.results_dir,
                    )

                    futures.append(f_chrn_df)

                for fut in futures:
                    res = fut.result()
                    future_chrn_df.append(res[1])
                    f_chrn_bench.append(res[0])

                benchmarks.extend(f_chrn_bench)

                for fut in future_chrn_df:
                    if individuals_filename in cfg.output_fns:
                        cfg.output_fns[individuals_filename].append(fut)
                    else:
                        cfg.output_fns[individuals_filename] = [fut]

                counter = counter + step

            # Sifting Job
            f_sifting = row[2]
            f_sifting = os.path.join(
                cfg.data_dir, "data", cfg.dataset, "sifting", f_sifting
            )

            f_sifted = executor.submit(
                sifting, inputfile=f_sifting, c=c_num, results_dir=cfg.results_dir
            )

            cfg.sifted_files.append(f_sifted)

        # merge task
        individuals_files = {}
        for key, val in cfg.output_fns.items():
            for data in val:
                if key in individuals_files:
                    if data[0] in individuals_files[key]:
                        individuals_files[key][data[0]].append(data[1])
                    else:
                        individuals_files[key][data[0]] = [data[1]]
                else:
                    individuals_files[key] = {data[0]: [data[1]]}

    cfg.sifted_files = [s.result() for s in cfg.sifted_files]
    benchmarks.extend([s[0] for s in cfg.sifted_files])
    print("Sifting completed")

    # Analyses jobs
    mutations = []
    frequencies = []

    if cfg.debug:
        individuals_files = list(individuals_files.items())[0:2]
    else:
        individuals_files = individuals_files.items()

    for i, inf in enumerate(individuals_files):
        for f_pop in cfg.populations:
            mutation_res = executor.submit(
                run_moverlap,
                input_dir=inf[1],
                siftfile=cfg.sifted_files[i][1],
                c=cfg.c_nums[i],
                data_dir=cfg.data_dir,
                results_dir=cfg.results_dir,
                columns=cfg.columns,
                pop=f_pop,
                debug=cfg.debug,
            )
            mutations.append(mutation_res)

            frequency_res = executor.submit(
                run_frequency,
                input_dir=inf[1],
                siftfile=cfg.sifted_files[i][1],
                c=cfg.c_nums[i],
                data_dir=cfg.data_dir,
                results_dir=cfg.results_dir,
                columns=cfg.columns,
                pop=f_pop,
                debug=cfg.debug,
            )
            frequencies.append(frequency_res)

    benchmarks.extend([m.result() for m in mutations])
    benchmarks.extend([freq.result() for freq in frequencies])
    return benchmarks


def run_proxy_workflow(cfg: Workflow, executor: Executor, store: Store) -> list[Bench]:
    benchmarks: list[Bench] = []

    with open(cfg.datafile, "r") as f:
        for row in csv.reader(f):
            base_file = row[0]
            threshold = int(cfg.fraction * int(row[1]))
            # To ensure we do not create too many individuals jobs
            cfg.ind_jobs = min(cfg.ind_jobs, threshold)
            step = threshold // cfg.ind_jobs
            rest = threshold % cfg.ind_jobs
            if rest != 0:
                raise ValueError(
                    f"For file {base_file}: required individuals jobs "
                    f"{cfg.ind_jobs} does not divide the number of rows "
                    f"{threshold}."
                )

            counter = 1

            output_files = []

            # Individuals Jobs
            f_individuals = base_file
            base_file_path = os.path.join(
                cfg.data_dir, "data", cfg.dataset, f_individuals
            )

            c_num = base_file[base_file.find("chr") + 3 :]
            c_num = c_num[0 : c_num.find(".")]
            cfg.c_nums.append(c_num)

            if cfg.debug:
                cfg.c_nums = [1, 2]
                threshold = 2

            individuals_filename = "chr%sn" % c_num

            while counter < threshold:
                stop = counter + step

                out_name = os.path.join("chr%sn-%s-%s.tar.gz" % (c_num, counter, stop))
                output_files.append(out_name)
                chrp_data_future: ProxyFuture[tuple] = store.future(
                    polling_interval=0.001
                )
                data_future: ProxyFuture[...] = store.future(polling_interval=0.001)

                chrom_bench = executor.submit(
                    processing_chrom_parts,
                    inputfile=base_file_path,
                    columnfile=os.path.join(
                        cfg.data_dir, "data", cfg.dataset, cfg.columns
                    ),
                    c=c_num,
                    counter=counter,
                    stop=stop,
                    total=threshold,
                    results_dir=cfg.results_dir,
                    chrp_data_future=chrp_data_future,
                    data_future=data_future,
                )

                ndir, chrp_elements = chrp_data_future.result()
                benchmarks.append(chrom_bench.result())

                print("Done obtaining chromosome parts")
                f_chrn_bench = []

                if cfg.debug:
                    chrp_elements = chrp_elements[0:5]

                futures = [
                    store.future(polling_interval=0.001)
                    for i in range(len(chrp_elements))
                ]
                # pf_chrn_df: ProxyFuture[(pd.DataFrame, str)] = store.future()

                for i, chrp in enumerate(chrp_elements):
                    f_chrn_bench.append(
                        executor.submit(
                            processing_2,
                            chrp_element=chrp,
                            data=data_future.proxy(),
                            ndir=ndir,
                            c=c_num,
                            results_dir=cfg.results_dir,
                            pfuture=futures[i],
                        )
                    )

                    name = f'chr{c_num}.{chrp["name"]}'
                    if individuals_filename not in cfg.output_fns:
                        cfg.output_fns[individuals_filename] = []
                    cfg.output_fns[individuals_filename].append(
                        (name, futures[i].proxy())
                    )

                benchmarks.extend(f_chrn_bench)

                counter = counter + step

            # Sifting Job
            f_sifting = row[2]
            f_sifting = os.path.join(
                cfg.data_dir, "data", cfg.dataset, "sifting", f_sifting
            )

            pf_sifting = store.future(polling_interval=0.001)
            f_sifted = executor.submit(
                sifting,
                inputfile=f_sifting,
                c=c_num,
                results_dir=cfg.results_dir,
                pfuture=pf_sifting,
            )

            cfg.sifted_files.append(pf_sifting.proxy())
            benchmarks.append(f_sifted)

        # merge task
        print("Merging results")
        individuals_files = {}
        for key, val in cfg.output_fns.items():
            for name, data_proxy in val:
                if key in individuals_files:
                    if name in individuals_files[key]:
                        individuals_files[key][name].append(data_proxy)
                    else:
                        individuals_files[key][name] = [data_proxy]
                else:
                    individuals_files[key] = {name: [data_proxy]}

    print("Sifting completed")

    # Analyses jobs
    mutations = []
    frequencies = []

    if cfg.debug:
        individuals_files = list(individuals_files.items())[0:2]
    else:
        individuals_files = individuals_files.items()

    for i, inf in enumerate(individuals_files):
        for f_pop in cfg.populations:
            mutation_res = executor.submit(
                run_moverlap,
                input_dir=inf[1],
                siftfile=cfg.sifted_files[i],
                c=cfg.c_nums[i],
                columns=cfg.columns,
                pop=f_pop,
                debug=cfg.debug,
            )
            mutations.append(mutation_res)

            frequency_res = executor.submit(
                run_frequency,
                input_dir=inf[1],
                siftfile=cfg.sifted_files[i],
                c=cfg.c_nums[i],
                columns=cfg.columns,
                pop=f_pop,
                debug=cfg.debug,
            )
            frequencies.append(frequency_res)

    benchmarks = [b.result() if isinstance(b, ComputeFuture) else b for b in benchmarks]
    benchmarks.extend([m.result() for m in mutations])
    print("Collected mutations")
    benchmarks.extend([freq.result() for freq in frequencies])
    print("Collected frequencies")
    return benchmarks


def main(argv: Sequence[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument(
        "--executor",
        choices=["globus-compute", "process-pool"],
        required=True,
    )
    parser.add_argument(
        "--proxystore",
        action="store_true",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
    )
    parser.add_argument(
        "--endpoint",
        help="Globus Compute Endpoint UUID",
    )
    parser.add_argument(
        "--fraction",
        default=1.0,
        type=float,
        help="Fraction of data in each file to use",
    )
    parser.add_argument(
        "--ind-jobs",
        default=16,
        type=int,
    )
    args = parser.parse_args()

    w = Workflow(
        data_dir=args.data_dir,
        results_dir=args.results_dir,
        ind_jobs=args.ind_jobs,
        debug=args.debug,
        fraction=args.fraction,
    )

    tic = perf_counter()

    store = Store("genome-store", FileConnector("/tmp/proxystore-cache"))

    if args.executor == "globus-compute":
        with Executor(endpoint_id=args.endpoint) as gce:
            if args.proxystore:
                benchmarks = run_proxy_workflow(w, executor=gce, store=store)
            else:
                benchmarks = run_gc_workflow(w, executor=gce)
    elif args.executor == "process-pool":
        with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as pool:
            if args.proxystore:
                benchmarks = run_proxy_workflow(w, executor=pool, store=store)
            else:
                benchmarks = run_gc_workflow(w, executor=pool)
    else:
        raise ValueError(f"Unsupported executor: {args.executor}")

    store.close()

    duration = perf_counter() - tic
    print(
        f"Workflow completed in {duration:.3f} (executor={args.executor}, "
        f"proxystore={args.proxystore}, debug={args.debug})"
    )

    name = args.executor
    if args.proxystore:
        name += "-proxystore"
    if args.debug:
        name += "-debug"

    df = process_benchmarks(benchmarks)
    df.to_pickle(f"{name}.pkl")

    if args.debug:
        plt.barh(
            y=df["thread_id"],
            width=df["duration"],
            left=df["norm_start"],
            color=df["colour"],
        )
        plt.savefig(f"{name}.pdf")


if __name__ == "__main__":
    raise SystemExit(main())
