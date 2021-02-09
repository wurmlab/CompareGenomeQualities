#!/bin/bash -l
# Here, -l is required to load conda environment. Without this we get:
#
#   CommandNotFoundError: Your shell has not been properly configured
#   to use 'conda activate'.
#
# Copyright 2021 Anurag Priyam, Queen Mary University of London.

usage() {
  cat <<EOF
CompareGenomeQualities

Usage:
compare-genome-qualities.sh -g 300000000 -b insecta_odb9 \
-1 illumina_R1.fq.gz -2 illumina_R2.fq.gz \
assembly_1.fa assembly_2.fa assembly_3.fa ...

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
EOF
  exit
}

msg() {
  echo >&2 -e "${1-}"
}

die() {
  local msg=$1
  local code=${2-1} # default exit status 1
  msg "$msg"
  exit "$code"
}

# Determine location of the script.
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)

# Default values for optional arguments.
num_cpus=1
output_dir=compare-genome-qualities-$(date +"%Y-%m-%d-%H%M%S")

# Parse arguments.
while :; do
  case "${1-}" in
  --rank-only)
    shift 1
    exec ${script_dir}/rank_assemblies.sh "$@"
    ;;
  -b | --busco-lineage)
    busco_lineage=${2}
    shift 2
    ;;
  -g | --genome-size)
    genome_size=${2}
    shift 2
    ;;
  -1 | --illumina-R1)
    illumina_r1=${2}
    shift 2
    ;;
  -2 | --illumina-R2)
    illumina_r2=${2}
    shift 2
    ;;
  -n | --num-cpus)
    num_cpus=${2}
    shift 2
    ;;
  -o | --output-dir)
    output_dir=${2}
    shift 2
    ;;
  -h | --help) usage ;;
  -?*) die "Unknown option: $1" ;;
  *) break ;;
  esac
done
assemblies=("$@")

# If we did not get all the required parameters, print usage and exit.
[[ ${#assemblies[@]} -eq 0 || -z $genome_size || -z $busco_lineage ]] && usage
[[ -n $illumina_r1 && -z $illumina_r2 ]] && usage
[[ -n $illumina_r2 && -z $illumina_r1 ]] && usage

# If Illumina reads were given, check the files exist, are not empty,
# and readable.
if [[ -n $illumina_r1 ]]; then
  [[ ! -s $illumina_r1 || ! -r $illumina_r1 ]] && die "${illumina_r1} is empty or not readable"
  [[ ! -s $illumina_r2 || ! -r $illumina_r2 ]] && die "${illumina_r2} is empty or not readable"
fi

# Check the given assembly files exist, are not empty and readable.
for file in "${assemblies[@]}"; do
  [[ ! -s ${file} || ! -r ${file} ]] && die "${file} is empty or not readable"
done

# Create output directory. -p is required because it may be desirable to write
# to the same output directory across multiple runs.
mkdir -p ${output_dir}

# To run BUSCO we need to download the lineage dataset first.
BUSCO_BASE_URL=https://busco-data.ezlab.org/v5/data/lineages
BUSCO_FILENAME=( $(curl -sL ${BUSCO_BASE_URL} | grep ${busco_lineage} | awk -F'>|<' '{print $3}') )
if [[ ${#BUSCO_FILENAME[@]} -eq 1 ]]; then
  echo "Downloading ${BUSCO_BASE_URL}/${BUSCO_FILENAME}"
  wget -c -P ${output_dir} ${BUSCO_BASE_URL}/${BUSCO_FILENAME} \
  && tar -C ${output_dir} -xvf ${output_dir}/${BUSCO_FILENAME} \
  || die "Error downloading BUSCO dataset."
else
  [[ ${#BUSCO_FILENAME[@]} -eq 0 ]] && die "Lineage ${busco_lineage} did not match any BUSCO dataset"
  [[ ${#BUSCO_FILENAME[@]} -gt 1 ]] && die "Lineage ${busco_lineage} matched multiple datasets: ${BUSCO_FILENAME[@]}"
fi

# Ensure the script exits if one of the subsequent steps fail.
set -eo pipefail

# Load parallel.
conda activate parallel

# Setup a function to run parallel with logging enabled. Note that we can't use
# aliases as they are not exported to subshell.
parallel() {
  command parallel --verbose --header : --results ${output_dir}/logs "$@" > /dev/null
}
export -f parallel

# Run QUAST.
(
  conda activate --stack quast
  parallel ${script_dir}/quast.sh -g ${genome_size} \
  -o ${output_dir}/'$(basename {})'/quast {} ::: quast.sh ${assemblies[@]}
)

# Run BUSCO.
(
  conda activate --stack busco
  parallel ${script_dir}/busco.sh -b ${output_dir}/${BUSCO_FILENAME%%.*} \
  -n ${num_cpus} -o ${output_dir}/'$(basename {})'/busco {} ::: busco.sh ${assemblies[@]}
)

# Read mapping analysis.
if [[ -n $illumina_r1 ]]; then
  # Map Illumina reads.
  (
    conda activate --stack bwa
    conda activate --stack samtools
    parallel ${script_dir}/map_illumina.sh -n ${num_cpus} \
    -o ${output_dir}/'$(basename {})'/illumina {} \
    ${illumina_r1} ${illumina_r2} ::: map_illumina.sh ${assemblies[@]}
  )

  # Calculate solid pairs and resolved length.
  (
    conda activate --stack ruby
    conda activate --stack samtools
    conda activate --stack bedtools
    conda activate --stack mosdepth
    parallel ${script_dir}/solid_pairs.sh \
    ${output_dir}/'$(basename {})'/illumina/bwa_mem_default.bam \
    ::: solid_pairs.sh ${assemblies[@]}
  )
fi

# Collect statistics.
${script_dir}/collect_stats.sh ${output_dir}

# Rank assemblies.
(
  conda activate --stack r
  ${script_dir}/rank_assemblies.sh ${output_dir}
)