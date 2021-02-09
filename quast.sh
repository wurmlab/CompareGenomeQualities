#!/bin/bash
# Copyright 2021 Anurag Priyam, Queen Mary University of London.

set -eo pipefail

usage() {
  cat <<EOF
quast.sh --genome-size value --output-dir value assembly.fa

CompareGenomeQualities wrapper script for QUAST.

--help              View this message
-o, --output-dir    Output directory. Required.
-g, --genome-size   Estimated genome size. Required.
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

while :; do
  case "${1-}" in
  -h | --help) usage ;;
  -g | --genome-size)
    genome_size=${2}
    shift 2
    ;;
  -o | --output-dir)
    output_dir=${2}
    shift 2
    ;;
  -?*) die "Unknown option: $1" ;;
  *) break ;;
  esac
done
assembly=$1

# If we did not get all the required parameters, print usage and exit.
[[ -z $genome_size || -z $output_dir || -z $assembly ]] && usage

# Check assembly file exists, is not-emptry, and readable.
if [[ ! -s ${assembly} || ! -r ${assembly} ]]; then
  die "${assembly} doesn't exist, or is empty, or not readable"
fi

# Run quast.
quast --fast -t 1 -e -m 0 --est-ref-size ${genome_size} \
  -o ${output_dir} ${assembly}
