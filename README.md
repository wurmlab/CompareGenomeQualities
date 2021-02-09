CompareGenomeQualities tool, described in “Priyam et al. (in review) Parameter exploration improves the accuracy of long-read genome assembly.”

CompareGenomeQualities ranks assemblies based on contiguity, completeness, and accuracy. Contiguity is assessed through NG50 metric. Completeness and accuracy are assessed in genic regions using BUSCO and genome-wide using short Illumina reads.

CompareGenomeQualities is distributed as a docker image containing all of its dependencies (GNU parallel, QUAST, BUSCO, bwa, samtools, bedtools, mosdepth, Ruby, R, numo-narray Ruby package, and ggplot2, ggpubr, rmarkdown, tidyr, and scales R packages). Docker is free to use and works on Windows, Linux, and Mac. You can install it from: https://www.docker.com/products/docker-desktop.

## Get tool and view help

```bash
docker run wurmlab/compare-genome-qualities --help

CompareGenomeQualities

Usage:
compare-genome-qualities.sh -g 300000000 -b insecta_odb9 -1 illumina_R1.fq.gz -2 illumina_R2.fq.gz assembly_1.fa assembly_2.fa assembly_3.fa ...

or,
compare-genome-qualities.sh --rank-only dir_containing_tsv_files

Options:
-b, --busco-lineage One of the 193 BUSCO v5 datasets listed here: https://busco-data.ezlab.org/v5/data/lineages.
                    Names can be a partial match e.g. insecta instead of insecta_odb10.2020-09-10.tar.gz.
                    Required unless --rank-only is specified.
-g, --genome-size   Expected or estimated genome size in base pairs. Required unless --rank-only is specified.
-1, --illumina-R1   Forward Illumina reads. Required if -2 is specified.
-2, --illumina-R2   Reverse Illumina reads. Required if -1 is specified.
-n, --num-cpus      Used for read mapping and BUSCO steps. Default: 1.
-o, --output-dir    Output directory. Default: `pwd`/compare-genome-qualities-yyyy-mm-dd-hhmmss.
                    Not applicable if --rank-only is specified.
--rank-only         Don’t compute metrics. Only rank assemblies based on tabular files in the given directory.
-h, --help          View this message
```

## Analysing the assemblies
To compare assemblies and select the best one, provide the path to the assemblies as input to the tool, along with the estimated genome size (for NG50 calculation), name of the BUSCO dataset to download and use (see https://busco-data.ezlab.org/v5/data/lineages), and the path to paired Illumina reads.


```bash
# --volume option makes the current directory available inside docker.
# input/ is a folder in the current directory containing assemblies
# and paired Illumina reads.
docker run --volume=`pwd`:/mnt:rw wurmlab/compare-genome-qualities \
--genome-size 450000000 \
--busco-lineage insecta \
--illumina-R1 input/illumina_R1.fq.gz \
--illumina-R2 input/illumina_R2.fq.gz \
input/assembly_1.fa input/assembly_2.fa
```

## Output
The tool generates intermediate output and final results to a compare-genome-qualities-timestamp/ folder in the current directory.

```bash
ls compare-genome-qualities-2020-12-23-150523/
assembly1/
assembly2/
.tab files
.Rmd
.pdf
```

## Including additional metrics

To include additional metrics, add a tabular file (see syntax below) to the output directory of CompareGenomeQualities and run:

```bash
docker run --volume=`pwd`:/mnt:rw wurmlab/compare-genome-qualities --rank-only compare-genome-qualities-2020-12-23-150523
```

The ranking mechanism can be used with a completely different set of metrics as well. The `--rank-only` option only expects a folder containing tabular files (see syntax below).

```bash
docker run -v $PWD:/mnt wurmlab/compare-genome-qualities --rank-assemblies my_metrics/
```

### Syntax of tabular files
One line per assembly. Two values on each line, separated by a tab.

```bash
assembly_id	metric_value
```
