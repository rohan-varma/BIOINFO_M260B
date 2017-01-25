import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from BIOINFO_M260B.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref
THRESHOLD = 7 # mismatches threshold

def mismatchify_pair(front, back, ref):
    """Takes a front/back read and matches it to a location in the genome.
       This algorithm slides the read through the entire genome, repeatedly saving
       the index and number of mismatches when it finds a smallest mismatching.
       Params:
        front - string, front read
        back: - string, back read
        ref: - string, large ref genome
       return:
        min_mismatches: min number of forward mismatches
        min_mismatch_loc: the location where the above was found in ref
        min_back_mismatches: the min number of backward mismatches
        min_back_loc: the location where the above was found in ref
        """
    min_mismatches, min_mismatch_loc = len(front) + 1, -1

    for i in range(len(ref) - len(front)):
        n_mismatches = sum([1 if front[j] != ref[i + j] else 0 for j in range(len(front))])
        if n_mismatches < min_mismatches:
            min_mismatches, min_mismatch_loc = n_mismatches, i

    if min_mismatches >= 10:
        return min_mismatches, min_mismatch_loc, 0, -1 # early termination

    min_back_mismatches, min_back_loc = len(back) + 1, -1

    for k in range(130,170):
        endpoint = min_mismatch_loc + k
        if endpoint + len(back) > len(ref): break # out of bounds
        n_mismatches = sum([1 if back[j] != ref[endpoint + j] else 0 for j in range(len(back))])
        if n_mismatches < min_back_mismatches:
            min_back_mismatches, min_back_loc = n_mismatches, endpoint

    return min_mismatches, min_mismatch_loc, min_back_mismatches, min_back_loc

def print_info(start, count, paired_end_reads):
    """Print details about time and how many reads have been aligned."""
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
    matched = 0
    # go through all reads
    for read_pair in paired_end_reads:
        count += 1
        if count % 50 == 0: print_info(start, count, paired_end_reads)
        front, reversed_back = read_pair[0], read_pair[1][::-1]
        # get mismatches
        min_mismatches, min_mismatch_loc, min_back_mismatches, min_back_loc = mismatchify_pair(front, reversed_back, ref)
        # check threshold. If met, add to list
        if (min_mismatches + min_back_mismatches < THRESHOLD):
            matched += 1
            read_alignment_locations = [min_mismatch_loc, min_back_loc]
            output_read_pair = [front, reversed_back]
            all_read_alignment_locations.append(read_alignment_locations)
            output_read_pairs.append(output_read_pair)

    # finish
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
