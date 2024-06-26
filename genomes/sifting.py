import logging
import os
import re
import sys
import threading
import time

import pandas as pd

from genomes.utils import Bench

logger = logging.getLogger(__name__)


def readfile(file):
    with open(file) as f:
        content = f.readlines()
    return content


def sifting(inputfile, c, results_dir, data_future=None, dask=False):
    tic = time.perf_counter()
    start = time.time()

    # unzipped = 'ALL.chr{}.vcf'.format(c)
    final = os.path.join(results_dir, f"sifted.SIFT.chr{c}.txt")

    rawdata = readfile(inputfile)

    logger.debug(f"= Taking columns from {inputfile}")
    logger.debug(f"== Filtering over {len(rawdata)} lines")

    # header =$(head - n 1000 $unzipped | grep "#" | wc - l)
    r1 = re.compile(".*(#).*")
    header = len(list(filter(r1.match, rawdata[:1000])))
    logger.debug(f"== Header found -> {header}")

    # grep -n "deleterious\|tolerated" $unzipped | grep "rs" > $siftfile
    ## This regex is super slow
    # r2 = re.compile('.*(deleterious|tolerated).*')
    data_temp = [
        f"{i}:{line}"
        for i, line in enumerate(rawdata)
        if "deleterious" in line or "tolerated" in line
    ]  # list(filter(r2.match, rawdata))
    # data_temp = []

    # init_size = len(rawdata)

    # for lineno,content in enumerate(rawdata[header:]):
    #     if r2.search(content):
    #         data_temp.append(str(lineno)+':'+content)
    #     logger.debug('{}/{}'.format(lineno, init_size), end='\r')

    # siftfile = 'SIFT.chr{}.vcf'.format(c)
    # with open(siftfile, 'w') as f:
    #     subprocess.run(["grep -n \"deleterious\|tolerated\" {}".format(inputfile)], shell=True, stdout=f)

    # data_temp = readfile(siftfile)

    r3 = re.compile(".*(rs).*")
    data = list(filter(r3.match, data_temp))

    logger.debug(f"== Starting processing {len(data)} lines")

    df = pd.DataFrame(columns=["line", "id", "ENSG_id", "SIFT", "phenotype"])

    for i, l in enumerate(data):
        # awk '{print $1}' $siftfile | awk -F ":" '{print $1-'$header'}' > $lines #.txt
        line = str(int(l.split("\t")[0].split(":")[0]) - int(header))
        # awk '{print $3}' $siftfile > $ids #.txt
        id = l.split("\t")[2]

        # awk '{print $8}' $siftfile > $info  # .txt
        # awk - F "|" '{print $5"\t"$17"\t"$18}' $info | sed 's/(/\t/g' | sed 's/)//g' > $sifts
        sifts = l.split("\t")[7].split("|")
        sifts = sifts[4] + " " + sifts[16] + " " + sifts[17]
        sifts = sifts.replace("(", " ").replace(")", "")

        # pr -m -t -s ' ' $lines $ids $sifts | gawk '{print $1,$2,$3,$5,$7}' > $final
        temp = (line + " " + id + " " + sifts).split(" ")

        if temp[3] == "" or temp[4] == "":
            df.loc[i] = [temp[0], temp[1], temp[2], None, None]
        elif temp[5] == "":
            df.loc[i] = [temp[0], temp[1], temp[2], temp[4], None]
        else:
            df.loc[i] = [temp[0], temp[1], temp[2], temp[4], temp[6]]

    if data_future is None and not dask:
        df.to_pickle(final)
    elif data_future is not None:
        data_future.set_result(df)
    # df.to_pickle(final)
    # if data_future is not None:
    #     data_future.set_result(final)

    logger.debug(
        f"= Line, id, ENSG id, SIFT, and phenotype printed to {final} in {time.perf_counter() - tic:0.2f} seconds.",
    )

    duration = time.perf_counter() - tic
    end = time.time()

    benchmark = Bench(threading.get_native_id(), "sifting", start, end, duration)

    if data_future is not None:
        return benchmark
    elif dask:
        return (benchmark, df.to_csv())
    return (benchmark, final)


if __name__ == "__main__":
    sifting(inputfile=sys.argv[1], c=sys.argv[2])
