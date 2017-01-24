import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from BIOINFO_M260B.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref

def process_read_pair(front, back, ref):
    min_mismatches, min_mismatch_location = len(front) + 1, -1
    for i in range(len(ref) - len(front)):
        n_mismatches = sum([1 if front[j] != ref[i + j] else 0 for j in range(len(front))])
        if n_mismatches < min_mismatches:
            min_mismatches, min_mismatch_location = n_mismatches, i

    min_tail_mismatches, min_tail_location = len(back) + 1, -1
    for k in range(130,170):
        endpoint = min_mismatch_location + k
        if endpoint + len(back) > len(ref): break # out of bounds
        n_mismatches = sum([1 if back[j] != ref[endpoint + j] else 0 for j in range(len(back))])
        if n_mismatches < min_tail_mismatches:
            min_tail_mismatches, min_tail_location = n_mismatches, endpoint

    return min_mismatches, min_mismatch_location, min_tail_mismatches, min_tail_location

def print_info(start, count, paired_end_reads):
    time_passed = (time.clock()-start)/60
    print('{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed))
    remaining_time = time_passed/count*(len(paired_end_reads)-count)
    print('Approximately {:.3} minutes remaining'.format(remaining_time))
    print count

def align_pairs_to_ref(paired_end_reads, ref):
    """
    Align paired end reads to the reference while minimizing errors.
    Params:
    paired_end_reads: list of (front, back) reads
    ref: A reference genome
    return:
    all_alignment_locations: list of alignment locations for each read
    output_read_pairs: list of readpairs oriented in their ideal orientation
    """
    all_read_alignment_locations, output_read_pairs = [], []
    count = 0
    start = time.clock()
    debug = False
    matched = 0
    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations, output_read_pairs = [], []
        found = False
        if count % 50 == 0: print_info(start, count, paired_end_reads)
        front, reversed_back = read_pair[0], read_pair[1][::-1]
        min_mismatches, min_mismatch_location, min_tail_mismatches, min_tail_location = process_read_pair(front, reversed_back, ref)

        if (min_mismatches + min_tail_mismatches < 10):
            matched += 1
            read_alignment_locations.append(min_mismatch_location)
            output_read_pair.append(front)
            read_alignment_locations.append(min_tail_location)
            output_read_pair.append(reversed_back)
            found = True
        if (found):
            all_read_alignment_locations.append(read_alignment_locations)
            output_read_pairs.append(output_read_pair)

    # after main loop
    print('Matched {} reads'.format(matched))
    return all_read_alignment_locations, output_read_pairs


def do_main():
    data_folder = 'hw1_W_2'
    input_folder = join('../data/', data_folder)
    f_base = '{}_chr_1'.format(data_folder)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
    start = time.clock()
    input_reads = read_reads(reads_fn)
    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    alignments, reads = align_pairs_to_ref(input_reads, reference)
    output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    output_fn = join(input_folder, 'aligned_{}.txt'.format(f_base))
    print("saving to " + output_fn)
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)
