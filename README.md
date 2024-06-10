# 1000Genomes Workflow

The 1000 genomes project provides a reference for human variation, having reconstructed the genomes of 2,504 individuals across 26 different populations to energize these approaches. This workflow identifies mutational overlaps using data from the 1000 genomes project in order to provide a null distribution for rigorous statistical evaluation of potential disease-related mutations. The workflow fetchs, parses, and analyzes data from the [1000 genomes Project](https://www.internationalgenome.org) (see ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/). It cross-matches the extracted data (which person has which mutations), with the mutation's sift score (how bad it is). Then it performs a few analyses, including plotting.

The figure below shows a branch of the workflow for the analysis of a single chromosome.

<p align="center">
  <img src="workflow.png?raw=true" width="600">
</p>

_Individuals_. This task fetches and parses the Phase3 data from the 1000 genomes project by chromosome. These files list all of Single nucleotide polymorphisms (SNPs) variants in that chromosome and which individuals have each one. SNPs are the most common type of genetic variation among people, and are the ones we consider in this work. An individual task creates output files for each individual of _rs numbers_ 3, where individuals have mutations on both alleles.

_Populations_. The 1000 genome project has 26 different populations from many different locations around the globe. A population task downloads a file per population selected. This workflow uses six super populations: African (AFR), Mixed American (AMR), East Asian (EAS), European (EUR), British from England and Scotland (GBR) and South Asian (SAS). The workflow also uses ALL population, which means that all individuals from the latest release are considered.

_Sifting_. A sifting task computes the SIFT scores of all of the SNPs variants, as computed by the Variant Effect Predictor (VEP). SIFT is a sequence homology-based tool that Sorts Intolerant From Tolerant amino acid substitutions, and predicts whether an amino acid substitution in a protein will have a phenotypic effect. For each chromosome the sifting task processes the corresponding VEP, and selects only the SNPs variants that has a SIFT score, recording in a file (per chromosome) the SIFT score and the SNPs variants ids, which are: (1) rs number, (2) ENSEMBL GEN ID, and (3) HGNC ID.

_Mutations_Overlap_. This task measures the overlap in mutations (also called SNPs variants) among pairs of individuals by population and by chromosome.

_Frequency_. This tasks measures the frequency of overlapping in mutations by selecting a number of random individuals, and selecting all SNPs variants without taking into account their SIFT scores.


This workflow is based on the application described in: https://github.com/rosafilgueira/Mutation_Sets

## Installation

Create a Python virtual environment or Conda environment and install the `genomes` package.
The package requires Python 3.8--3.11, but may work with newer versions.

```bash
python -m venv venv
. venv/bin/activate
pip install -e .
```

## Prepare the Data

```bash
./prepare-data.sh
```

## Run the Application
```
./daxgen.py -D 20130502 -f data.csv -i 1
```
This workflow assumes that all input data listed in the `data.csv` file is available in the `data/20130502` folder by default (but you can change that behavior with the `-D`).

### Workflow parallelism
You can control how many `individuals` jobs **per chromosome** will get created with the parameter `-i IND_JOBS, --individuals-jobs IND_JOBS`, by default it's set to `1`. If the value provided is larger than the total number of rows in the data file for that chromosome, then it will be set to the number of rows so that each job will process one row (_Warning_: this will extremely inefficient and will create a large number of jobs, about `250,000`).

In addition, it is required that `IND_JOBS` **divides the number of rows** for each chromosome, in this case `250,000`.

# References

Dreher, Matthieu, and Tom Peterka. _Decaf: Decoupled dataflows for in situ high-performance workflows._ No. ANL/MCS-TM-371. Argonne National Lab.(ANL), Argonne, IL (United States), 2017. https://www.mcs.anl.gov/~tpeterka/papers/2017/dreher-anl17-report.pdf
