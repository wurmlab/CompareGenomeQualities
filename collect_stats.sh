#!/bin/bash
# Copyright 2021 Anurag Priyam, Queen Mary University of London.

set -eo pipefail

usage() {
  cat <<EOF
CompareGenomeQualities script for collecting assembly metrics.

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

# Parse arguments.
while :; do
  case "${1-}" in
  -h | --help) usage ;;
  -?*) die "Unknown option: $1" ;;
  *) break ;;
  esac
done
directory=$1

# If we did not get all the required parameters, print usage and exit.
[[ -z $directory ]] && usage

# We need to cd to the given directory.
cd ${directory}

# Output files.
ng50_tsv=NG50.tsv
busco_tsv=busco.tsv
solid_tsv=solid_pairs.tsv
resolved_tsv=resolved_length.tsv

exists() {
  compgen -G "$1" > /dev/null
}

# NG50
if exists '*/quast/report.txt'; then
  echo -e "ID\tNG50" > ${ng50_tsv}
  grep NG50 */quast/report.txt \
  | cat | awk -F'/|[ \t]+' -v OFS='\t' '{print $1,$4}' >> ${ng50_tsv}
fi

# BUSCO score
if exists '*/busco/busco_output/short_summary*.txt'; then
  echo -e "ID\tBUSCO" > ${busco_tsv}
  grep 'C:' */busco/busco_output/short_summary*.txt \
  | cat | awk -F'[]/%:[ \t]+' -v OFS='\t' '{print $1,$6}' >> ${busco_tsv}
fi

if exists '*/illumina/solid_pairs.txt'; then
  # Solid pairs
  echo -e "ID\tSolid pairs" > ${solid_tsv}
  grep 'Solidly mapped' */illumina/solid_pairs.txt \
  | cat | awk -F'[/:)( \t]+' -v OFS='\t' '{print $1,$9}' >> ${solid_tsv}

  # Resolved length
  echo -e "ID\tResolved length" > ${resolved_tsv}
  grep 'Resolved length' */illumina/solid_pairs.txt \
  | cat | awk -F'[/: \t]+' -v OFS='\t' '{print $1,$6}' >> ${resolved_tsv}
fi