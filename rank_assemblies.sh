#!/bin/bash
# Copyright 2021 Anurag Priyam, Queen Mary University of London.

usage() {
  cat <<EOF
CompareGenomeQualities script for ranking assemblies.

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

# Determine location of the script.
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)

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

output_dir=${directory}

cp -r ${script_dir}/rank_assemblies.Rmd ${output_dir}
R -e "rmarkdown::render('${output_dir}/rank_assemblies.Rmd', 'all')"