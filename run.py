from __future__ import annotations

import argparse
import csv
import logging
import multiprocessing
import os
import sys
import threading
import time
from concurrent.futures import Future as ComputeFuture
from concurrent.futures import ProcessPoolExecutor
from typing import Sequence

import pandas as pd
from globus_compute_sdk import Executor
from proxystore.connectors.file import FileConnector
from proxystore.store import Store

from genomes.frequency import run_frequency
from genomes.individuals import processing_chrom_parts
from genomes.mutation_overlap import run_moverlap
from genomes.sifting import sifting
from genomes.utils import Bench

logger = logging.getLogger(__name__)


def process_benchmarks(benchmarks):
    df = pd.DataFrame(benchmarks)
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
    fraction: float

    def __init__(
        self,
        data_dir: str,
        results_dir: str,
        ind_jobs: int = 1,
        fraction: float = 1.0,
    ) -> None:
        self.data_dir = data_dir
        self.results_dir = results_dir
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
    individuals_files = {}

    with open(cfg.datafile) as f:
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
                    f"{threshold}.",
                )

            # Individuals Jobs
            f_individuals = base_file
            base_file_path = os.path.join(
                cfg.data_dir, "data", cfg.dataset, f_individuals,
            )

            c_num = base_file[base_file.find("chr") + 3 :]
            c_num = c_num[0 : c_num.find(".")]
            cfg.c_nums.append(c_num)

            individuals_filename = "chr%sn" % c_num

            logger.info(f"Processing individual: {individuals_filename}")
            future_chrn_df = []

            counter = 1
            while counter < threshold:
                stop = counter + step

                logger.info(f"Submitting parts: {counter} - {stop}")
                task_future = executor.submit(
                    processing_chrom_parts,
                    inputfile=base_file_path,
                    columnfile=os.path.join(
                        cfg.data_dir, "data", cfg.dataset, cfg.columns,
                    ),
                    c=c_num,
                    counter=counter,
                    stop=stop,
                    total=threshold,
                    results_dir=cfg.results_dir,
                )
                future_chrn_df.append(task_future)

                counter = counter + step

            output_fns = {}
            for fut in future_chrn_df:
                if individuals_filename not in output_fns:
                    output_fns[individuals_filename] = []
                output_fns[individuals_filename].append(fut)

            logger.info("Completed individuals job")

            merge_start = time.time()
            # merge task
            logger.info("Merging results")
            for key, futures in output_fns.items():
                for future in futures:
                    bench, results = future.result()
                    benchmarks.append(bench)
                    for name, df in results:
                        if key in individuals_files:
                            if name in individuals_files[key]:
                                individuals_files[key][name].append(df)
                            else:
                                individuals_files[key][name] = [df]
                        else:
                            individuals_files[key] = {name: [df]}
            merge_end = time.time()
            benchmarks.append(
                Bench(
                    threading.get_native_id(),
                    "merge",
                    merge_start,
                    merge_end,
                    merge_end - merge_start,
                ),
            )

            # Sifting Job
            f_sifting = row[2]
            f_sifting = os.path.join(
                cfg.data_dir, "data", cfg.dataset, "sifting", f_sifting,
            )

            f_sifted = executor.submit(
                sifting, inputfile=f_sifting, c=c_num, results_dir=cfg.results_dir,
            )

            cfg.sifted_files.append(f_sifted)

            logger.info("Completed sifting")

    cfg.sifted_files = [s.result() for s in cfg.sifted_files]
    benchmarks.extend([s[0] for s in cfg.sifted_files])
    logger.info("Sifting completed")

    # Analyses jobs
    mutations = []
    frequencies = []

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
            )
            frequencies.append(frequency_res)

    benchmarks.extend([m.result() for m in mutations])
    benchmarks.extend([freq.result() for freq in frequencies])
    return benchmarks


def run_proxy_workflow(cfg: Workflow, executor: Executor, store: Store) -> list[Bench]:
    benchmarks: list[Bench] = []
    individuals_files = {}

    with open(cfg.datafile) as f:
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
                    f"{threshold}.",
                )

            # Individuals Jobs
            f_individuals = base_file
            base_file_path = os.path.join(
                cfg.data_dir, "data", cfg.dataset, f_individuals,
            )

            c_num = base_file[base_file.find("chr") + 3 :]
            c_num = c_num[0 : c_num.find(".")]
            cfg.c_nums.append(c_num)

            individuals_filename = "chr%sn" % c_num

            logger.info(f"Processing individual: {individuals_filename}")
            individual_task_futures = []
            individual_data_futures = []

            counter = 1
            while counter < threshold:
                stop = counter + step

                logger.info(f"Submitting parts: {counter} - {stop}")
                data_future = store.future(polling_interval=0.001)
                task_future = executor.submit(
                    processing_chrom_parts,
                    inputfile=base_file_path,
                    columnfile=os.path.join(
                        cfg.data_dir, "data", cfg.dataset, cfg.columns,
                    ),
                    c=c_num,
                    counter=counter,
                    stop=stop,
                    total=threshold,
                    results_dir=cfg.results_dir,
                    data_future=data_future,
                )
                individual_task_futures.append(task_future)
                individual_data_futures.append(data_future)

                counter = counter + step

            for data_future in individual_data_futures:
                data = data_future.proxy()
                if individuals_filename not in cfg.output_fns:
                    cfg.output_fns[individuals_filename] = []
                cfg.output_fns[individuals_filename].append(data)

            benchmarks.extend(individual_task_futures)

            logger.info("Completed individuals job")

            # Sifting Job
            f_sifting = row[2]
            f_sifting = os.path.join(
                cfg.data_dir, "data", cfg.dataset, "sifting", f_sifting,
            )

            data_future = store.future(polling_interval=0.001)
            task_future = executor.submit(
                sifting,
                inputfile=f_sifting,
                c=c_num,
                results_dir=cfg.results_dir,
                data_future=data_future,
            )

            cfg.sifted_files.append(data_future.proxy())
            benchmarks.append(task_future)

            logger.info("Completed sifting")

        merge_start = time.time()
        # merge task
        logger.info("Merging results")
        for key, proxy in cfg.output_fns.items():
            for val in proxy:
                for name, df in val:
                    if key in individuals_files:
                        if name in individuals_files[key]:
                            individuals_files[key][name].append(df)
                        else:
                            individuals_files[key][name] = [df]
                    else:
                        individuals_files[key] = {name: [df]}
        merge_end = time.time()
        benchmarks.append(
            Bench(
                threading.get_native_id(),
                "merge",
                merge_start,
                merge_end,
                merge_end - merge_start,
            ),
        )

    # Analyses jobs
    mutations = []
    frequencies = []

    individuals_files = individuals_files.items()

    for i, inf in enumerate(individuals_files):
        for f_pop in cfg.populations:
            mutation_res = executor.submit(
                run_moverlap,
                input_dir=inf[1],
                siftfile=cfg.sifted_files[i],
                c=cfg.c_nums[i],
                data_dir=cfg.data_dir,
                results_dir=cfg.results_dir,
                columns=cfg.columns,
                pop=f_pop,
            )
            mutations.append(mutation_res)

            frequency_res = executor.submit(
                run_frequency,
                input_dir=inf[1],
                siftfile=cfg.sifted_files[i],
                c=cfg.c_nums[i],
                data_dir=cfg.data_dir,
                results_dir=cfg.results_dir,
                columns=cfg.columns,
                pop=f_pop,
            )
            frequencies.append(frequency_res)

    benchmarks = [b.result() if isinstance(b, ComputeFuture) else b for b in benchmarks]
    benchmarks.extend([m.result() for m in mutations])
    benchmarks.extend([freq.result() for freq in frequencies])
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
    parser.add_argument(
        "--workers",
        type=int,
    )
    args = parser.parse_args()

    log_fmt = "[%(asctime)s.%(msecs)03d] %(levelname)-5s (%(name)s) :: %(message)s"
    logging.basicConfig(
        format=log_fmt,
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )

    w = Workflow(
        data_dir=args.data_dir,
        results_dir=args.results_dir,
        ind_jobs=args.ind_jobs,
        fraction=args.fraction,
    )

    tic = time.perf_counter()

    store = Store("genome-store", FileConnector("/tmp/proxystore-cache"))

    if args.executor == "globus-compute":
        with Executor(endpoint_id=args.endpoint, batch_size=8) as gce:
            if args.proxystore:
                benchmarks = run_proxy_workflow(w, executor=gce, store=store)
            else:
                benchmarks = run_gc_workflow(w, executor=gce)
    elif args.executor == "process-pool":
        workers = (
            args.workers if args.workers is not None else multiprocessing.cpu_count()
        )
        with ProcessPoolExecutor(max_workers=workers) as pool:
            if args.proxystore:
                benchmarks = run_proxy_workflow(w, executor=pool, store=store)
            else:
                benchmarks = run_gc_workflow(w, executor=pool)
    else:
        raise ValueError(f"Unsupported executor: {args.executor}")

    store.close()

    duration = time.perf_counter() - tic
    logger.info(
        f"Workflow completed in {duration:.3f} (executor={args.executor}, "
        f"proxystore={args.proxystore})",
    )

    name = args.executor
    if args.proxystore:
        name += "-proxystore"

    df = process_benchmarks(benchmarks)
    df["name"] = name
    filepath = f"run-{name}.csv"
    df.to_csv(filepath, index=True)
    logger.info(f"Saved task times to {filepath}")


if __name__ == "__main__":
    raise SystemExit(main())
