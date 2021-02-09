#!/bin/bash
# Copyright 2021 Anurag Priyam, Queen Mary University of London.

set -eo pipefail

usage() {
  cat <<EOF
CompareGenomeQualities wrapper script for BUSCO.

-b, --busco-dataset Path to BUSCO dataset to use. Required.
-o, --output-dir    Output directory. Required.
-n, --num-cpus      Number of CPUs to use for
                    all except BLAST stage.
--help              View this message
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

# Default values for optional arguments.
num_cpus=1

# Parse arguments.
while :; do
  case "${1-}" in
  -h | --help) usage ;;
  -b | --busco-dataset)
    busco_dataset=${2}
    shift 2
    ;;
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
assembly=$1

# If we did not get all the required arguments, print usage and exit.
[[ -z $busco_dataset || -z $output_dir || -z $assembly ]] && usage

# Create output directory
mkdir -p ${output_dir}

# BUSCO writes to current working directory so we expand input paths,
# cd to output directory and work from there.
busco_dataset=$(readlink -f ${busco_dataset})
assembly=$(readlink -f ${assembly})
cd ${output_dir}

# Run BUSCO.
busco --offline -m geno -c ${num_cpus} -l ${busco_dataset} -i ${assembly} -o busco_output