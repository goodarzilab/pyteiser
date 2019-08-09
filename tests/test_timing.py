# see how it's implemented here: https://www.kennethreitz.org/essays/repository-structure-and-python

import os
import sys
import argparse
import numpy as np
import timeit
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
            current_sequence = structures.s_sequence(len(sequence_string))
            current_sequence.from_sequence(sequence_string)

            time_create_object = timeit.timeit(lambda: structures.s_sequence(len(sequence_string)), number=100)
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

# def func():
#     var1 = 'aaa'
#     var2 = 'aab'
#     def closure():
#         return var1 == var2
#     t1 = timeit.timeit(closure, number = 10**4)
#     print(t1)


def main():
    args = handler()
    # func()
    time_reading_fasta(args.rna_fastafile)
    # time_compressing_sequences(args.rna_fastafile)









if __name__ == "__main__":
    main()

