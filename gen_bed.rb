#!/usr/bin/env ruby
# Copyright 2019 Anurag Priyam, Queen Mary U London

# The input BAM file.
input_bam = ARGV.shift

# Contig length threshold.
min_length = ARGV.shift

# Exclude ends by this much.
ignore_ends = ARGV.shift

# Print usage and exit with status code 1 unless all three parameters
# have been provided.
unless input_bam && min_length && ignore_ends
  fail "Usage: gen_bed.rb input_bam min_contig_length ignore_ends"
  exit!
end

# Convert min_length and ignore_ends to integer.
min_length = min_length.to_i
ignore_ends = ignore_ends.to_i

# Read in reference sequence sizes from the BAM header.
IO.foreach("| samtools view -H #{input_bam} | grep '^@SQ'") do |line|
  ref, len = line.split.values_at(1, 2)
  ref = ref[3..-1]
  len = len[3..-1].to_i
  if len >= min_length
    sart = ignore_ends
    stop = len - ignore_ends
    puts [ref, ignore_ends, len - ignore_ends].join("\t")
  end
end
