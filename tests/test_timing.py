# see how it's implemented here: https://www.kennethreitz.org/essays/repository-structure-and-python

import os
import sys
import argparse
import numpy as np
import timeit
import numba
import bitarray
import random
sys.path.insert(0, os.path.abspath('..'))

import pyteiser.glob_var as glob_var
import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.matchmaker as matchmaker



def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seedfile", type=str)
    parser.add_argument("--rna_fastafile", type=str)
    parser.add_argument("--outfile", type=str)


    parser.set_defaults(
        rna_fastafile = '/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/narrow_down_transcripts_list/Gencode_v28_GTEx_expressed_transcripts_fasta/utr_3_fasta/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.txt',
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
    )

    args = parser.parse_args()

    return args

# passing callables to timeit: use labmda
# read here: https://stackoverflow.com/questions/31550832/timeit-timeit-variable-importing-in-python



def time_reading_fasta(fasta_file):
    tr_dict_loc = {}
    seqs_order = []
    with open(fasta_file, 'r') as f:
        split_string = f.read().split('>')
        for entry in split_string:
            if entry == '':
                continue
            seq_start = entry.find('\n')
            annotation = entry[:seq_start]
            sequence_string = entry[seq_start + 1:].replace('\n', '')
            current_sequence = structures.w_sequence(len(sequence_string))
            current_sequence.from_sequence(sequence_string)

            time_create_object = timeit.timeit(lambda: structures.w_sequence(len(sequence_string)), number=100)
            time_fill_object = timeit.timeit(lambda: current_sequence.from_sequence(sequence_string), number=100)
            time_compress_object = timeit.timeit(lambda: current_sequence.compress(), number=100)
            time_compress_named_object = timeit.timeit(lambda: IO.compress_named_sequences({annotation: current_sequence}, [annotation]), number=100)

            print("Create object: %.5f" % time_create_object)
            print("Fill object: %.5f" % time_fill_object)
            print("Compress object: %.5f" % time_compress_object)
            print("Compress named object: %.5f" % time_compress_named_object)
            print()


            # curr_timing = timeit.timeit('current_sequence.from_sequence(sequence_string)',
            #                             'from __main__ import current_sequence, sequence_string')
            # print(curr_timing)

            #
            # tr_dict_loc[annotation] = current_sequence
            # seqs_order.append(annotation)

    return tr_dict_loc, seqs_order


def time_compressing_sequences(fasta_file):
    sequences_dict, seqs_order = IO.read_fasta(fasta_file)

    for i in range(len(seqs_order)):
        print(seqs_order[i])


@numba.jit(cache=True, nopython=True, nogil=True)
def do_iterate(x):
    a = 0
    for i in x:
        a = 1


def time_iterating():
    pr_list = [0] * 100000
    pr_numpy = np.array(pr_list)
    time_just_list = timeit.timeit(lambda: do_iterate(pr_list), number=10)
    time_numpy_list = timeit.timeit(lambda: do_iterate(pr_numpy), number=10)

    print("Iterating through 100k long list 10 times takes: ", time_just_list)
    print("Iterating through 100k numpy array 10 times takes: ", time_numpy_list)

    # iterating through numpy array takes 50 times faster than iterating through a list


def time_compressing_profile():
    TOY_ARRAY_LENGTH = 10000
    NUMBER_OF_ONES = 2000

    ones_indices_list = random.sample(range(TOY_ARRAY_LENGTH), NUMBER_OF_ONES)
    ones_indices_set = set(ones_indices_list)
    toy_array = [1 if x in ones_indices_set else 0 for x in range(TOY_ARRAY_LENGTH)]

    toy_bool_array = np.asarray(toy_array, dtype=bool)

    toy_bitarray = bitarray.bitarray(toy_array)
    toy_packbits_array = np.packbits(toy_bool_array)

    time_bitarray = timeit.timeit(lambda: bitarray.bitarray(toy_array), number=10000)
    time_packbits = timeit.timeit(lambda: np.packbits(toy_bool_array), number=10000)

    print("Bitarray conversion takes", time_bitarray)
    print("Numpy conversion takes", time_packbits)

    # Numpy packbits is ~13 times faster than bitarray


def main():
    args = handler()

    # time_reading_fasta(args.rna_fastafile)
    # time_compressing_sequences(args.rna_fastafile)
    # time_iterating()
    # time_compressing_profile()




if __name__ == "__main__":
    main()

