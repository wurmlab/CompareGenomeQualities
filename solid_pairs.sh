#!/bin/bash
# Copyright 2021 Anurag Priyam, Queen Mary University of London.

set -eo pipefail

usage() {
  cat <<EOF
CompareGenomeQualities solid pairs metric.

Usage:
solid_pairs.sh bwa_mem_default.bam

Options:
--window-size       Window size for calculations. Default: 1 bp
--output-dir        Output directory. Default: same directory
                    as input BAM file
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

# Assume window size of 1bp; can be overriden (-w).
window_size=1
while :; do
  case "${1-}" in
  -h | --help) usage ;;
  -w | --window-size)
    window_size=${$2}
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
genome_bam=${1}

# If we did not get all the required parameters, print usage and exit.
[[ -z $genome_bam ]] && usage

# Set output directory if none specified.
[[ -z $output_dir ]] && output_dir=$(dirname ${genome_bam})

# Determine location of the script.
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)

# Output prefix.
prefix=${output_dir}/${window_size}bp

# Intermediate and final output files.
lowconf_bed=${prefix}.regions.lowconf.bed
highconf_bed=${prefix}.regions.highconf.bed
highconf_bam=${prefix}.regions.highconf.bam
solid_pairs_txt=${prefix}.solid_pairs.txt

# Get mean depth in a window.
mosdepth -xn -b ${window_size} ${prefix} ${genome_bam}

# Decompress mosdepth's output as we are going to read it twice.
# -k causes the compressed file to be retained (default behaviour
# is to delete the source file).
gunzip -kf ${prefix}.regions.bed.gz

# Calculate and output mean, median and modal read depth of the windows.
${script_dir}/mean_mode_median.rb ${prefix}.regions.bed > ${solid_pairs_txt}

# Eliminate regions coverage less than 5x coverage or higher than twice the
# median coverage.
median=$(awk -F', ' '{ print $6 }' ${solid_pairs_txt})
${script_dir}/filter_windows.rb ${prefix}.regions.bed \
5 $(( ${median} * 2 )) > ${highconf_bed} 2> ${lowconf_bed}

echo -e "\nAfter filtering ...\n" >> ${solid_pairs_txt}

# Output resolved length.
num_highconf_regions=$(wc -l ${highconf_bed} | awk '{print $1}')
echo "Resolved length: $(( num_highconf_regions * window_size ))" \
>> ${solid_pairs_txt}

# Output mean, median, mode of resolved regions. If the mean equals
# median and the mode, it is a good sign.
${script_dir}/mean_mode_median.rb ${highconf_bed} >> ${solid_pairs_txt}

# -sorted switch reduces the memory usage considerably. This requires the files
# to be sorted: mosdepth's BED output is sorted by default and we sort BAM as a
# standard processing step.
bedtools intersect -sorted -a ${genome_bam} -b ${highconf_bed} > ${highconf_bam}

# Count and output solidly mapped Illumina read pairs.
num_solid_pairs=$(ruby solid_pairs.rb ${highconf_bam})
num_reads=$(samtools view -F 2304 ${genome_bam} | wc -l | awk '{print $1}')
perc_solid_pairs=$(perl -e "printf('%.2f', 2 * ${num_solid_pairs} * 100 / ${num_reads})")
echo "Solidly mapped read pairs: ${num_solid_pairs} (${perc_solid_pairs})" >> ${solid_pairs_txt}

# Cleanup
# 1. This files takes a lot of space
rm ${prefix}.regions.bed
# 2. Human friendly output file name
ln -sf $(basename ${solid_pairs_txt}) ${output_dir}/solid_pairs.txt
