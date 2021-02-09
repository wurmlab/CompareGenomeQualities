#!/usr/bin/env ruby
# Copyright 2019 Anurag Priyam, Queen Mary U London
#
# Script to calculate mean, median, and modal coverage values from BED-like
# with coverage in last column. 

require 'numo/narray'

# The input BAM file (compulsory).
input_bed = ARGV.shift

# Check the file exists and is not empty.
if !File.exist?(input_bed) || File.zero?(input_bed)
  puts 'File does not exist, or it is empty.'
  exit!
end

# Determine number of lines in the BED file.
num_windows = IO.foreach(input_bed).count

# Read depth values from the bed file.
depth = Numo::UInt16.zeros(num_windows)
i = 0
IO.foreach(input_bed) do |line|
  depth[i] = line.split.last.to_i
  i += 1
end

mean = (depth.sum / num_windows).round
mode = depth.bincount.max_index
median = depth.median

sd = 0
depth.each do |d|
  sd += (mean - d) ** 2
end
sd = Math.sqrt(sd / num_windows).round

puts "Mean, median, mode, sd: #{mean}, #{median}, #{mode}, #{sd}"
