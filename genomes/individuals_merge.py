import logging
import os
import shutil
import sys
import tarfile
import tempfile
import time

logger = logging.getLogger(__name__)


def compress(archive, input_dir):
    with tarfile.open(archive, "w:gz") as f:
        f.add(input_dir, arcname="")


def extract_all(archive, output_dir):
    with tarfile.open(archive, "r:*") as f:
        f.extractall(os.path.basename(output_dir))
        flist = f.getnames()
        if flist[0] == "":
            flist = flist[1:]
        return flist


def readfile(filename):
    with open(filename) as f:
        content = f.readlines()
    return content


def writefile(filename, content):
    with open(filename, "w") as f:
        f.writelines(content)


def merging(c, tar_files):
    logger.debug(f"= Merging chromosome {c}...")
    tic = time.perf_counter()

    merged_dir = f"merged_chr{c}"
    os.makedirs(merged_dir, exist_ok=True)

    data = {}

    for tar in tar_files:
        tic_iter = time.perf_counter()
        with tempfile.TemporaryDirectory(dir=os.curdir) as temp_dir:
            for filename in extract_all(tar, temp_dir):
                content = readfile(os.path.join(temp_dir, filename))
                if filename in data:
                    data[filename] += content
                else:
                    data[filename] = content

        logger.debug(
            f"Merged {tar} in {time.perf_counter() - tic_iter:0.2f} sec",
        )

    for filename, content in data.items():
        writefile(os.path.join(merged_dir, filename), content)

    outputfile = f"chr{c}n.tar.gz"
    logger.debug(f"== Done. Zipping {len(data)} files into {outputfile}.")

    compress(outputfile, merged_dir)

    # Cleaning temporary files
    try:
        shutil.rmtree(merged_dir)
    except OSError as e:
        logger.error("Error: %s : %s" % (merged_dir, e.strerror))

    logger.debug(
        f"= Chromosome {c} merged in {time.perf_counter() - tic:0.2f} seconds.",
    )


if __name__ == "__main__":
    print(f"Host = {os.uname()[1]}")
    print(f"CPUs = {os.sched_getaffinity(0)}")
    merging(c=sys.argv[1], tar_files=sys.argv[2:])
