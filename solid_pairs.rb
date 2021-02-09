#!/usr/bin/env ruby
# Copyright 2019 Anurag Priyam, Queen Mary U London
#
# Prints the number of solid read pairs in the given BAM file. By
# solid we mean reads that mapped end-to-end and in a proper pair.

# The input BAM file.
input_bam = ARGV.shift

# Check the file exists and is not empty.
if !File.exist?(input_bam) || File.zero?(input_bam)
  puts 'File does not exist, or it is empty.'
  exit!
end

# Returns true if the read mapped in a proper pair and the read is not clipped.
def read_passed?(flag, cigar)
  proper_pair?(flag) && not_clipped?(cigar)
end

# Returns true if the read mapped in a proper pair. This does not check if the
# mate is present in the BAM file. So it will return true even if the read was
# orphaned by a preprocessing step, such as bedtools intersect step of the cmg
# pipeline. This is okay. The read is eliminated later on as not a solid pair.
def proper_pair?(flag)
  flag.to_i.to_s(2)[-2..-1] == '11'
end

# Returns true if the read is not clipped.
def not_clipped?(cigar)
  cigar !~ /(\d+)S/
end

# The first time a read is encountered it's id and pass/fail status is stored in
# a lookup table. The second time a read is encountered, this lookup table is
# consulted and the fate of the pair decided. If a read is not encountered a
# second time, it must have been orphaned upstream.
lookup_table = {}

# Number of solidly mapping read pairs.
num_solid_pairs = 0

# -F 2304 removes secondary and supplementary alignments. bwa-mem does not
# output secondary alignments by default, but other mappers may. We filter
# these because we only want to consider one representative alignment of
# each read - whichever the mapper chose as representative.
IO.foreach("| samtools view -F 2304 #{input_bam}") do |line|
  # We need fragment id, flag, and cigar string for each read.
  id, flag, cigar = line.split.values_at(0, 1, 5)

  # If we have not seen this read before, add its pass or fail status to the
  # lookup table. But if we are seeing a read id for the second time, it must be
  # the mate and we can now decide about the pair as a whole. If the read mapped
  # as a proper pair but was orphaned by a prior step, we will never encounter
  # its mate and it will not count as solid pair.
  if !lookup_table.include?(id)
    lookup_table[id] = read_passed?(flag, cigar)
  else
    num_solid_pairs += 1 if lookup_table.delete(id) && read_passed?(flag, cigar)
  end
end

print num_solid_pairs