#!/usr/bin/env python3

import time

tic = time.perf_counter()
import numpy as np
from random import sample
import os.path
import matplotlib
import pandas as pd
import tarfile
import threading

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
import collections
from collections import Counter
from io import StringIO

from genomes.utils import Bench

from proxystore.proxy import extract, Proxy

SIFT = "NO-SIFT"


class ReadData:
    def __init__(
        self,
        input_dict: dict[str, dict[str, list]],
        pop_dir: str,
        data_dir: str,
        chrom: str,
        debug: bool = False,
    ):
        self.input_dict = input_dict
        self.pop_dir = pop_dir
        self.data_dir = data_dir
        self.chrom = chrom
        self.debug = debug

    def read_names(self, POP):
        print("reading inidviduals")
        tic = time.perf_counter()
        namefile = self.pop_dir + POP
        f = open(namefile, "r")
        text = f.read()
        f.close()
        text = text.split()
        all_ids = text[0:]
        file = self.data_dir + "columns.txt"
        f = open(file, "r")
        text = f.read()
        f.close()
        genome_ids = text.split()

        ids = list(set(all_ids) & set(genome_ids))
        print("time: %s" % (time.perf_counter() - tic))
        return ids

    def read_rs_numbers(self, siftfile):
        print("reading in rs with sift scores below %s" % SIFT)
        ## NB This file is in the format of:
        ## line number, rs number, ENSG number, SIFT, Phenotype
        tic = time.perf_counter()
        rs_numbers = []
        variations = {}
        map_variations = {}
        all_variations = []

        try:
            if os.path.exists(siftfile):
                sift_file = pd.read_pickle(siftfile)
            elif isinstance(siftfile, str):
                sift_file = pd.read_csv(StringIO(siftfile))
            else:
                sift_file = siftfile
        except TypeError:
            sift_file = siftfile

        for idx, item in sift_file.iterrows():
            if len(item) > 2:
                rs_numbers.append(item.iloc[2])
                map_variations[item.iloc[2]] = item.iloc[3]
                variations[item.iloc[1]] = item.iloc[3]

        print("time: %s" % (time.perf_counter() - tic))
        print("read_rs_numbers completed")
        return rs_numbers, map_variations

    def read_individuals(self, ids, rs_numbers):
        print("reading in individual mutation files")
        tic = time.perf_counter()
        mutation_index_array = []

        if self.debug:
            ids = list(set(n.split(".")[1] for n in self.input_dict.keys()))

        for name in ids:
            # filename = self.data_dir + chrom + 'n/' + chrom + '.' + name

            chrom_dir = f"{self.chrom}n"
            fn = f"{self.chrom}.{name}"

            merged_text = []
            for i in range(len(self.input_dict[fn])):
                df = self.input_dict[fn][i]

                if isinstance(df, str):
                    if df[0] != ",":
                        df = pd.read_pickle(df)
                    else:
                        df = pd.read_csv(StringIO(df))

                text = df["data"].to_list()
                text = [e for t in text for e in t.split()]
                merged_text.extend(text)

            sifted_mutations = list(set(rs_numbers).intersection(merged_text))
            mutation_index_array.append(sifted_mutations)

        print("time: %s" % (time.perf_counter() - tic))
        return mutation_index_array


class Results:
    def __init__(self, n_runs: int, n_indiv: int):
        self.n_runs = n_runs
        self.n_indiv = n_indiv

    def overlap_ind(self, ids, mutation_index_array):
        n_p = len(mutation_index_array)
        print(
            "calculating the number overlapings mutations between %s individuals selected randomly"
            % n_p
        )
        tic = time.perf_counter()
        list_p = np.linspace(0, n_p - 1, n_p).astype(int)
        mutation_overlap = []
        random_indiv = []
        for run in range(self.n_runs):
            randomized_list = sample(list(list_p), n_p)
            result = Counter()
            r_ids = []
            for pq in range(self.n_indiv):
                if 2 * pq >= len(randomized_list):
                    break
                b_multiset = collections.Counter(
                    mutation_index_array[randomized_list[2 * pq]]
                )
                print("time, inidividual: %s" % ids[randomized_list[2 * pq]])
                r_ids.append(ids[randomized_list[2 * pq]])
                result = result + b_multiset
            random_indiv.append(r_ids)
            mutation_overlap.append(result)
        print("time: %s" % (time.perf_counter() - tic))
        return mutation_overlap, random_indiv

    def histogram_overlap(self, mutation_overlap):
        print("calculating the frequency/historgram of overlapings mutations")
        tic = time.perf_counter()
        histogram_overlap = []
        for run in range(self.n_runs):
            final_counts = [count for item, count in mutation_overlap[run].items()]
            histogram_overlap.append(collections.Counter(final_counts))
        print("time: %s" % (time.perf_counter() - tic))
        return histogram_overlap


class PlotData:
    def __init__(self, n_runs: int):
        self.n_runs = n_runs

    def plot_histogram_overlap(self, POP, histogram_overlap, outputFile):
        print("ploting Histogram mutation overlap to %s" % outputFile)
        tic = time.perf_counter()
        for run in range(self.n_runs):
            output = outputFile + str(run) + ".png"
            final_counts = [count for item, count in histogram_overlap[run].items()]
            N = len(final_counts)
            x = range(N)
            width = 1 / 1.5
            bar1 = plt.bar(x, final_counts, width, color="grey")
            plt.ylabel("Mutations")
            plt.xlabel("Individuals")
            plt.xticks(np.arange(1, N + 1))
            plt.savefig(output)
            plt.close()
        print("time: %s" % (time.perf_counter() - tic))


class WriteData:
    def __init__(self, n_runs: int, n_indiv: int):
        self.n_runs = n_runs
        self.n_indiv = n_indiv

    def write_histogram_overlap(self, histogram_overlapfile, histogram_overlap):
        print(
            "writing Frequency historgram of mutations overlapping to %s"
            % histogram_overlapfile
        )
        tic = time.perf_counter()
        for run in range(self.n_runs):
            overlapfile = histogram_overlapfile + str(run) + ".txt"
            f = open(overlapfile, "w")
            f.write("Number Individuals - Number Mutations  \n")
            for i in range(1, self.n_indiv + 1):
                if i in histogram_overlap[run]:
                    f.write(str(i) + "-" + str(histogram_overlap[run][i]) + "\n")
                else:
                    f.write(str(i) + "-" + str(0) + "\n")
            f.close()

        # for key, count in histogram_overlap[run].items() :
        # f.write(str(key) + '-' + str(count) + '\n')
        print("time: %s" % (time.perf_counter() - tic))

    def write_mutation_overlap(self, mutation_overlapfile, mutation_overlap):
        print("writing Mutations overlapping to %s" % mutation_overlapfile)
        tic = time.perf_counter()
        for run in range(self.n_runs):
            overlapfile = mutation_overlapfile + str(run) + ".txt"
            f = open(overlapfile, "w")
            f.write("Mutation Index- Number Overlapings \n")
            for key, count in mutation_overlap[run].items():
                f.write(key + "-" + str(count) + "\n")
            f.close()
        print("time: %s" % (time.perf_counter() - tic))

    def write_random_indiv(self, randomindiv_file, random_indiv):
        tic = time.perf_counter()
        for run in range(self.n_runs):
            randomfile = randomindiv_file + str(run) + ".txt"
            f = open(randomfile, "w")
            print("writing Random individuals to %s" % randomfile)
            f.write("Individuals \n")
            for item in random_indiv[run]:
                f.write("%s\n" % item)
            f.close()
        print("time: %s" % (time.perf_counter() - tic))

    def write_mutation_index_array(
        self, mutation_index_array_file, mutation_index_array
    ):
        print("writing Mutation index array to %s" % mutation_index_array_file)
        tic = time.perf_counter()
        f = open(mutation_index_array_file, "w")
        for item in mutation_index_array:
            f.write("%s\n" % item)
        f.close()
        print("time: %s" % (time.perf_counter() - tic))

    def write_map_variations(self, map_variations_file, map_variations):
        print("writing map_variations to %s" % map_variations_file)
        tic = time.perf_counter()
        f = open(map_variations_file, "w")
        for key, count in map_variations.items():
            f.write(key + "\t" + str(count) + "\n")
        f.close()
        print("time: %s" % (time.perf_counter() - tic))


############################################################
def run_frequency(
    input_dir,
    siftfile,
    c,
    pop,
    data_dir,
    results_dir,
    columns=None,
    debug: bool = False,
):
    tic = time.perf_counter()
    start = time.time()

    POP = pop
    n_runs = 1000
    n_indiv = 52

    chrom = "chr" + str(c)
    base_data_dir = data_dir
    data_dir = os.path.join(base_data_dir, "data/20130502/")
    pop_dir = os.path.join(base_data_dir, "data/populations/")
    outdata_dir = os.path.join(
        results_dir,
        "chr{0}-{1}-freq/output_no_sift/".format(str(c), str(POP)),
    )
    plot_dir = os.path.join(
        results_dir,
        "chr{0}-{1}-freq/plots_no_sift/".format(str(c), str(POP)),
    )

    if not os.path.exists(outdata_dir):
        os.makedirs(outdata_dir, exist_ok=True)
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir, exist_ok=True)

    OutputFormat = ".png"

    font = {"family": "serif", "size": 14}
    plt.rc("font", **font)

    rd = ReadData(
        input_dict=input_dir,
        pop_dir=pop_dir,
        data_dir=data_dir,
        chrom=chrom,
        debug=debug,
    )
    res = Results(n_runs=n_runs, n_indiv=n_indiv)
    wr = WriteData(n_runs=n_runs, n_indiv=n_indiv)
    pd = PlotData(n_runs=n_runs)

    histogram_overlapfile = (
        outdata_dir
        + "Histogram_mutation_overlap_chr"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + "_"
    )
    mutation_overlapfile = (
        outdata_dir
        + "Mutation_overlap_chr"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + "_"
    )
    mutation_index_array_file = (
        outdata_dir
        + "mutation_index_array"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + ".txt"
    )
    histogram_overlap_plot = (
        plot_dir + "Frequency_mutations" + str(c) + "_s" + str(SIFT) + "_" + POP
    )
    map_variations_file = (
        outdata_dir + "map_variations" + str(c) + "_s" + str(SIFT) + "_" + POP + ".txt"
    )

    randomindiv_file = (
        outdata_dir + "random_indiv" + str(c) + "_s" + str(SIFT) + "_" + POP + "_"
    )

    ids = rd.read_names(POP)
    n_pairs = len(ids) / 2

    if isinstance(siftfile, Proxy):
        siftfile = extract(siftfile)
    rs_numbers, map_variations = rd.read_rs_numbers(siftfile)
    mutation_index_array = rd.read_individuals(ids, rs_numbers)

    wr.write_map_variations(map_variations_file, map_variations)
    wr.write_mutation_index_array(mutation_index_array_file, mutation_index_array)

    mutation_overlap, random_indiv = res.overlap_ind(ids, mutation_index_array)
    histogram_overlap = res.histogram_overlap(mutation_overlap)

    wr.write_mutation_overlap(mutation_overlapfile, mutation_overlap)
    wr.write_histogram_overlap(histogram_overlapfile, histogram_overlap)
    wr.write_random_indiv(randomindiv_file, random_indiv)

    pd.plot_histogram_overlap(POP, histogram_overlap, histogram_overlap_plot)

    outfn = os.path.join(results_dir, "chr%s-%s-freq.tar.gz" % (c, POP))
    # gen final output
    tar = tarfile.open(outfn, "w:gz")
    tar.add(outdata_dir)
    tar.add(plot_dir)
    tar.close()

    duration = time.perf_counter() - tic
    end = time.time()
    return Bench(threading.get_native_id(), "frequency", start, end, duration)


if __name__ == "__main__":
    c_help = "type a chromosome 1-22"
    pop_help = "type a population 0-6; 0:ALL, 1:EUR, 2:EAS, 3:AFR, 4:AMR, 5:SAS, 6:GBR"
    description = "Process mutation sets (-c and -POP are required)."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-c", type=int, help=c_help)
    parser.add_argument("-pop", help=pop_help)
    args = parser.parse_args()
    c = args.c

    n_runs = 1000
    n_indiv = 52

    siftfile = "./sifted.SIFT.chr" + str(c) + ".txt"

    # untar input data
    import tarfile

    chrom = "chr" + str(c)
    input_dir = f"./{chrom}n"
    tar = tarfile.open(chrom + "n.tar.gz")
    tar.extractall(path=input)
    tar.close()

    run_frequency(input_dir=input_dir, siftfile=siftfile, c=c, pop=args.pop)
