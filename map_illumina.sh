#!/bin/bash
# Copyright 2021 Anurag Priyam, Queen Mary University of London.

set -eo pipefail

usage() {
  cat <<EOF
map_illumina.sh assembly.fa illumina_R1.fq.gz illumina_R2.fq.gz

CompareGenomeQualities wrapper script for mapping paired Illumina
reads using bwa-mem.

--help              View this message
-o, --output-dir    Output directory. Required.
-n, --num-cpus      Number of CPUs to use minus one. Required.
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

# ENV
TMPDIR=${TMPDIR:-/tmp}

# Parse arguments.
while :; do
  case "${1-}" in
  -h | --help) usage ;;
  -o | --output-dir)
    output_dir=${2}
    shift 2
    ;;
  -n | --num-cpus)
    num_cpus=${2}
    shift 2
    ;;
  -?*) die "Unknown option: $1" ;;
  *) break ;;
  esac
done
assembly=${1}
illumina_r1=${2}
illumina_r2=${3}

# If we did not get all the required parameters, print usage and exit.
[[ -z $output_dir || -z $num_cpus || -z $assembly || -z $illumina_r1 || -z $illumina_r2 ]] && usage

# Create output directory.
mkdir -p ${output_dir}

# Link assembly file to the output directory.
assembly_link=${output_dir}/$(basename ${assembly})
ln -s $(readlink -f ${assembly}) ${assembly_link}

# Index the assembly for read mapping.
bwa index ${assembly_link}

# Name of the BAM file we will create.
alignments_file="${output_dir}/bwa_mem_default.bam"

# Map reads and covert the resulting SAM output to sorted BAM format.
bwa mem -t ${num_cpus} ${assembly_link} ${illumina_r1} ${illumina_r2} \
2> ${output_dir}/bwa_mem.log | samtools sort -O BAM -@ ${num_cpus} \
-T ${TMPDIR} > ${alignments_file}

# Index the sorted BAM file.
samtools index -@ ${num_cpus} ${alignments_file}