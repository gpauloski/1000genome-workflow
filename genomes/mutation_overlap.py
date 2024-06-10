import argparse
import itertools
import logging
import os
import os.path
import tarfile
import threading
import time
from io import StringIO
from random import sample

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import pyplot
from proxystore.proxy import Proxy, extract

from genomes.utils import Bench

logger = logging.getLogger(__name__)

SIFT = "NO-SIFT"
n_runs = 1


class ReadData:
    def __init__(
        self,
        input_dict: dict[str, dict[str, list]],
        pop_dir: str,
        data_dir: str,
        chrom: str,
        debug=False,
    ):
        self.input_dict = input_dict
        self.pop_dir = pop_dir
        self.data_dir = data_dir
        self.chrom = chrom
        self.debug = debug

    def read_names(self, POP):
        logger.debug("reading individuals")
        tic = time.perf_counter()
        namefile = os.path.join(self.pop_dir, POP)
        f = open(namefile)
        text = f.read()
        f.close()
        text = text.split()
        all_ids = text[0:]
        file = self.data_dir + "columns.txt"
        f = open(file)
        text = f.read()
        f.close()
        genome_ids = text.split()

        ids = list(set(all_ids) & set(genome_ids))

        logger.debug("time: %s" % (time.perf_counter() - tic))
        return ids

    def read_rs_numbers(self, siftfile):
        logger.debug("reading in rs with sift scores below %s" % SIFT)
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

        for i, item in sift_file.iterrows():
            if len(item) > 3:
                rs_numbers.append(item.iloc[2])
                map_variations[item.iloc[2]] = item.iloc[3]

        # seen = set()
        # [x for x in all_variations if x not in seen and not seen.add(x)]
        # all_variations = seen

        logger.debug("time: %s" % (time.perf_counter() - tic))
        return rs_numbers, map_variations

    def read_individuals(self, ids, rs_numbers):
        logger.debug("reading in individual mutation files")
        tic = time.perf_counter()
        mutation_index_array = []
        total_mutations = {}
        total_mutations_list = []

        if self.debug:
            ids = list(set(n.split(".")[1] for n in self.input_dict.keys()))

        for name in ids:
            filename = self.data_dir + self.chrom + "n/" + self.chrom + "." + name

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
            total_mutations[name] = len(sifted_mutations)
            total_mutations_list.append(len(sifted_mutations))
            # logger.debug len(list(seen)), len(seen)

        logger.debug(
            "mutation index array for %s : %s" % (ids[0], mutation_index_array[0]),
        )
        logger.debug(
            "total_len_mutations for %s : %s" % (ids[0], total_mutations[ids[0]]),
        )
        logger.debug("total_mutations_list is %s " % total_mutations_list)
        logger.debug("time: %s" % (time.perf_counter() - tic))
        return mutation_index_array, total_mutations, total_mutations_list

    def read_pairs_overlap(self, indpairsfile):
        logger.debug("reading in individual crossover mutations")
        tic = time.perf_counter()
        pairs_overlap = np.loadtxt(indpairsfile, unpack=True)
        pairs_overlap = np.transpose(pairs_overlap)

        logger.debug("time: %s" % (time.perf_counter() - tic))
        return pairs_overlap


class Results:
    def __init__(self, n_runs: int, debug: bool = False) -> None:
        self.n_runs = n_runs
        self.debug = debug

    def group_indivuals(self, total_mutations_list):
        logger.debug("histograms mutations_individuals groups by 26")
        tic = time.perf_counter()

        n_group = 26
        if self.debug:
            n_group = 2

        random_mutations_list = []
        for run in range(self.n_runs):
            random_mutations_list.append(sample(total_mutations_list, n_group))
        logger.debug("time: %s" % (time.perf_counter() - tic))
        return random_mutations_list

    def pair_individuals(self, mutation_index_array):
        logger.debug("cross matching mutations in individuals")
        tic = time.perf_counter()

        n_p = len(mutation_index_array)
        n_pairs = int(round(n_p / 2))
        list_p = np.linspace(0, n_p - 1, n_p).astype(int)
        pairs_overlap = np.zeros((self.n_runs, n_pairs))
        for run in range(self.n_runs):
            randomized_list = sample(list(list_p), n_p)
            for pq in range(n_pairs):
                array1 = mutation_index_array[randomized_list[2 * pq]]
                array2 = mutation_index_array[randomized_list[2 * pq]]
                pair_array = set(array1) & set(array2)
                pairs_overlap[run][pq] = len(pair_array)

        logger.debug("time: %s" % (time.perf_counter() - tic))
        return pairs_overlap

    def total_pair_individuals(self, mutation_index_array):
        logger.debug("cross matching mutations total individuals")
        tic = time.perf_counter()
        n_p = len(mutation_index_array)
        total_pairs_overlap = np.zeros((n_p, n_p))
        simetric_overlap = np.zeros((n_p, n_p))
        for run in range(n_p):
            array1 = mutation_index_array[run]
            start = run + 1
            for pq in range(start, n_p):
                array2 = mutation_index_array[pq]
                pairs_array = set(array1) & set(array2)
                total_pairs_overlap[run][pq] = len(pairs_array)
                simetric_overlap[run][pq] = len(pairs_array)
                simetric_overlap[pq][run] = len(pairs_array)

        logger.debug("time: %s" % (time.perf_counter() - tic))
        return total_pairs_overlap, simetric_overlap

    def half_pair_individuals(self, mutation_index_array):
        logger.debug("cross matching mutations in individuals - half with half")
        tic = time.perf_counter()
        n_p = len(mutation_index_array)
        n_pairs = int(round(n_p / 2))
        pairs_overlap = np.zeros((n_pairs, n_pairs))
        for run in range(n_pairs):
            array1 = mutation_index_array[run]
            index = 0
            for pq in range(n_pairs + 1, n_p):
                array2 = mutation_index_array[pq]
                pairs_array = set(array1) & set(array2)
                pairs_overlap[run][index] = len(pairs_array)

        logger.debug("time: %s" % (time.perf_counter() - tic))
        return pairs_overlap

    def gene_pairs(self, mutation_index_array):
        logger.debug("cross matching pairs of variations")

        tic = time.perf_counter()
        n_p = len(mutation_index_array)
        gene_pair_list = {}
        for pp in range(n_p):
            pairs = itertools.combinations(mutation_index_array[pp], 2)
            for pair in pairs:
                key = str(pair)
                if key not in gene_pair_list:
                    gene_pair_list[key] = 1
                else:
                    gene_pair_list[key] += 1

        logger.debug("time: %s" % (time.perf_counter() - tic))

        return gene_pair_list


class PlotData:
    def __init__(self, c) -> None:
        self.c = c

    def individual_overlap(self, POP, pairs_overlap, outputFile):
        logger.debug(
            "plotting cross matched number of individuals:%s " % len(pairs_overlap),
        )
        tic = time.perf_counter()

        pairs_overlap = np.array(pairs_overlap)

        min_p = np.min(pairs_overlap)
        max_p = np.max(pairs_overlap)
        nbins = int(max_p) + 1
        n_runs = len(pairs_overlap)

        nbins = int(np.max(pairs_overlap))
        bin_centres = np.linspace(0, nbins, nbins)
        bin_edges = np.linspace(-0.5, nbins + 0.5, nbins + 1)

        fig = plt.figure(frameon=False, figsize=(10, 9))
        ax = fig.add_subplot(111)
        # ax.set_xlim([0,35])
        # ax.set_xlim([0,1000])
        hists = []
        max_h = 0
        for run in range(n_runs):
            h, edges = np.histogram(pairs_overlap[run], bins=bin_edges)
            ax.plot(bin_centres, h, alpha=0.5)
            if len(h) > 0:
                max_h = max(max_h, max(h))

        plt.xlabel("Number of overlapping gene mutations", fontsize=24)
        plt.ylabel(r"frequency", fontsize=28)
        text1 = (
            "population "
            + POP
            + "\n"
            + "chromosome "
            + str(self.c)
            + "\n"
            + "SIFT < "
            + str(SIFT)
            + "\n"
            + str(n_runs)
            + " runs"
        )
        plt.text(
            0.95,
            0.95,
            text1,
            fontsize=24,
            verticalalignment="top",
            horizontalalignment="right",
            transform=ax.transAxes,
        )
        plt.savefig(outputFile)
        plt.close()
        logger.debug("time: %s" % (time.perf_counter() - tic))

    def total_colormap_overlap(self, POP, total_pairs_overlap, outputFile):
        logger.debug(
            "plotting colormap number of individuals: %s" % len(total_pairs_overlap),
        )
        tic = time.perf_counter()
        fig = plt.figure()
        cmap = mpl.colors.ListedColormap(["blue", "black", "red", "green", "pink"])
        img = pyplot.imshow(
            total_pairs_overlap, interpolation="nearest", cmap=cmap, origin="lower",
        )
        pyplot.colorbar(img, cmap=cmap)

        # cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', ['blue','black','red'], 256)
        # img2 = pyplot.imshow(total_pairs_overlap,interpolation='nearest', cmap = cmap2, origin='lower')
        # pyplot.colorbar(img2,cmap=cmap2)

        # ax = fig.add_subplot(111)
        # ax.set_title('colorMap Cross Matched')
        # plt.imshow(total_pairs_overlap)
        # ax.set_aspect('equal')
        # cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
        # cax.get_xaxis().set_visible(False)
        # cax.get_yaxis().set_visible(False)
        # cax.patch.set_alpha(0)
        # cax.set_frame_on(False)
        # plt.colorbar(orientation='vertical')
        # plt.show()

        # x=list(range(0,len(total_pairs_overlap)))
        # y=list(range(0,len(total_pairs_overlap)))
        # x, y = np.meshgrid(x, y)
        # total_pairs_overlap=np.array(total_pairs_overlap)
        # plt.pcolormesh(x, y, total_pairs_overlap)
        # plt.colorbar() #need a colorbar to show the intensity scale
        # plt.show() #boom

        # sns.heatmap(total_pairs_overlap, vmax= max_p, yticklabels=50,  xticklabels=50)
        plt.savefig(outputFile)
        plt.close()
        logger.debug("time: %s" % (time.perf_counter() - tic))


class WriteData:
    def __init__(self, n_runs) -> None:
        self.n_runs = n_runs

    def write_pair_individuals(self, indpairsfile, pairs_overlap):
        logger.debug("writing pairs overlapping mutations to %s" % indpairsfile)
        tic = time.perf_counter()
        np.savetxt(indpairsfile, pairs_overlap, fmt="%i")
        logger.debug("time: %s" % (time.perf_counter() - tic))

    def write_gene_pairs(self, genepairsfile, gene_pair_list):
        logger.debug("writing gene pair list to %s" % genepairsfile)
        tic = time.perf_counter()
        f = open(genepairsfile, "w")
        for key, count in gene_pair_list.items():
            f.write(key + "\t" + str(count) + "\n")
        f.close()
        logger.debug("time: %s" % (time.perf_counter() - tic))

    def write_total_indiv(self, total_mutations_filename, total_mutations):
        logger.debug(
            "writing total mutations list per individual to %s"
            % total_mutations_filename,
        )
        tic = time.perf_counter()
        f = open(total_mutations_filename, "w")
        for key, count in total_mutations.items():
            f.write(key + "\t" + str(count) + "\n")
        f.close()
        logger.debug("time: %s" % (time.perf_counter() - tic))

    def write_random_mutations_list(
        self, random_mutations_filename, random_mutations_list,
    ):
        logger.debug(
            "writing a list of 26 random individuals with the number mutations per indiv %s"
            % random_mutations_filename,
        )
        for run in range(self.n_runs):
            filename = random_mutations_filename + "_run_" + str(run) + ".txt"
            f = open(filename, "w")
            f.writelines(["%s\n" % item for item in random_mutations_list[run]])

    def write_mutation_index_array(
        self, mutation_index_array_file, mutation_index_array,
    ):
        logger.debug("writing mutation array  to %s" % mutation_index_array_file)
        f = open(mutation_index_array_file, "w")
        for item in mutation_index_array:
            f.write("%s\n" % item)
        f.close()

    def write_map_variations(self, map_variations_file, map_variations):
        logger.debug("writing map_variations to %s" % map_variations_file)
        tic = time.perf_counter()
        f = open(map_variations_file, "w")
        for key, count in map_variations.items():
            f.write(key + "\t" + str(count) + "\n")
        f.close()


def run_moverlap(
    input_dir,
    siftfile,
    c,
    pop,
    data_dir,
    results_dir,
    columns=None,
    debug=False,
    proxy=None,
):
    tic = time.perf_counter()
    start = time.time()
    POP = pop
    chrom = "chr" + str(c)
    base_data_dir = data_dir
    data_dir = os.path.join(base_data_dir, "data/20130502/")
    pop_dir = os.path.join(base_data_dir, "data/populations/")
    outdata_dir = os.path.join(
        results_dir,
        f"chr{c!s}-{pop!s}/output_no_sift/",
    )
    plots_dir = os.path.join(
        results_dir,
        f"chr{c!s}-{pop!s}/plots_no_sift/",
    )

    if not os.path.exists(outdata_dir):
        os.makedirs(outdata_dir, exist_ok=True)
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir, exist_ok=True)

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
    res = Results(n_runs=n_runs, debug=debug)
    wr = WriteData(n_runs=n_runs)
    pd = PlotData(c=c)

    half_indpairsfile = (
        outdata_dir
        + "individual_half_pairs_overlap_chr"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + ".txt"
    )
    total_indpairsfile = (
        outdata_dir
        + "total_individual_pairs_overlap_chr"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + ".txt"
    )
    genepairsfile = (
        outdata_dir
        + "gene_pairs_count_chr"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + ".txt"
    )
    random_indpairsfile = (
        outdata_dir
        + "100_individual_overlap_chr"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + ".txt"
    )

    colormap = (
        plots_dir
        + "colormap_distribution_c"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + OutputFormat
    )
    half_overlap = (
        plots_dir
        + "half_distribution_c"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + OutputFormat
    )
    total_overlap = (
        plots_dir
        + "total_distribution_c"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + OutputFormat
    )
    random_overlap = (
        plots_dir
        + "100_distribution_c"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + OutputFormat
    )

    total_mutations_filename = (
        outdata_dir
        + "total_mutations_individual"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
        + ".txt"
    )
    random_mutations_filename = (
        outdata_dir
        + "random_mutations_individual"
        + str(c)
        + "_s"
        + str(SIFT)
        + "_"
        + POP
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

    map_variations_file = (
        outdata_dir + "map_variations" + str(c) + "_s" + str(SIFT) + "_" + POP + ".txt"
    )

    ids = rd.read_names(POP)
    n_pairs = len(ids) / 2

    if isinstance(siftfile, Proxy):
        siftfile = extract(siftfile)
    rs_numbers, map_variations = rd.read_rs_numbers(siftfile)
    mutation_index_array, total_mutations, total_mutations_list = rd.read_individuals(
        ids, rs_numbers,
    )
    wr.write_total_indiv(total_mutations_filename, total_mutations)
    wr.write_map_variations(map_variations_file, map_variations)

    # cross-correlations mutations overlapping
    half_pairs_overlap = res.half_pair_individuals(mutation_index_array)
    total_pairs_overlap, simetric_overlap = res.total_pair_individuals(
        mutation_index_array,
    )
    random_pairs_overlap = res.pair_individuals(mutation_index_array)

    wr.write_mutation_index_array(mutation_index_array_file, mutation_index_array)
    wr.write_pair_individuals(half_indpairsfile, half_pairs_overlap)
    wr.write_pair_individuals(total_indpairsfile, total_pairs_overlap)
    wr.write_pair_individuals(random_indpairsfile, random_pairs_overlap)

    pd.individual_overlap(POP, half_pairs_overlap, half_overlap)
    pd.individual_overlap(POP, simetric_overlap, total_overlap)
    pd.individual_overlap(POP, random_pairs_overlap, random_overlap)
    pd.total_colormap_overlap(POP, total_pairs_overlap, colormap)

    # list of frecuency of mutations in 26 individuals
    random_mutations_list = res.group_indivuals(total_mutations_list)
    wr.write_random_mutations_list(random_mutations_filename, random_mutations_list)

    # gen overlapping
    gene_pair_list = res.gene_pairs(mutation_index_array)
    wr.write_gene_pairs(genepairsfile, gene_pair_list)

    outfn = os.path.join(results_dir, "chr%s-%s.tar.gz" % (c, POP))
    # gen final output
    tar = tarfile.open(outfn, "w:gz")
    tar.add(outdata_dir)
    tar.add(plots_dir)
    tar.close()

    duration = time.perf_counter() - tic
    end = time.time()
    return Bench(threading.get_native_id(), "mutation_overlap", start, end, duration)


if __name__ == "__main__":
    c_help = "type a chromosome 1-22"
    pop_help = "type a population 0-6; 0:ALL, 1:EUR, 2:EAS, 3:AFR, 4:AMR, 5:SAS, 6:GBR"
    description = "Process mutation sets (-c and -POP are required)."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-c", type=int, help=c_help)
    parser.add_argument("-pop", help=pop_help)
    args = parser.parse_args()
    c = args.c

    n_runs = 1
    siftfile = "./sifted.SIFT.chr" + str(c) + ".txt"
    chrom = "chr" + str(c)

    # untar input data
    input_dir = f"./{chrom}n"
    import tarfile

    tar = tarfile.open(chrom + "n.tar.gz")
    tar.extractall(path="./" + chrom + "n")
    tar.close()

    run_moverlap(input_dir=input_dir, siftfile=siftfile, c=c, pop=args.pop)
