#!/usr/bin/env ruby
# Copyright 2019 Anurag Priyam, Queen Mary U London
#
# Script to filter windows by depth. There is both, a low-pass and a high-
# pass filter. Low-pass threshold is based on a priori expectation that
# windows with less number of supporting reads are likely to have higher
# sequencing errors. High-pass threshold is derived from the data: it is
# twice the median coverage.
#
# Low-coverage regions are ideally decided by mapping the reads used for
# assembly. By using Illumina read mappings I am assuming that the same
# regions will show up as low-coverage for both datasets. Illumina reads
# are ideally suited for high-pass filter. The key thing here is to not
# filter the BAM based on mapping quality.

require 'numo/narray'

# The input BAM file (compulsory).
input_bed = ARGV.shift
low_pass = Integer(ARGV.shift)
high_pass = Integer(ARGV.shift)

# Check the file exists and is not empty.
if !File.exist?(input_bed) || File.zero?(input_bed)
  puts 'File does not exist, or it is empty.'
  exit!
end

IO.foreach(input_bed) do |line|
  depth = line.split.last.to_i
  if depth >= low_pass && depth < high_pass
    puts line 
  else
    $stderr.puts line
  end
end
