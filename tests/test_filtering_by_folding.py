import argparse
import numpy as np
import timeit
import hashlib
import struct
import numba

import os
import sys

# to make sure relative import works
# for a detailed explanation, see test_matchmaker.py

current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)

import pyteiser.glob_var as glob_var
import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.matchmaker as matchmaker



def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seedfile", type=str)
    parser.add_argument("--rna_fastafile", type=str)
    parser.add_argument("--profiles_full_file", type=str)
    parser.add_argument("--profiles_filtered_file", help="", type=str)
    parser.add_argument("--indices_mode", help="compression in the index mode", type=bool)

    parser.set_defaults(
        rna_fastafile = '/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/narrow_down_transcripts_list/Gencode_v28_GTEx_expressed_transcripts_fasta/utr_3_fasta/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.txt',
        rna_bin_file='/Users/student/Documents/hani/temp/temp_bins/test_bin.bin',
        profiles_full_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique.bin',
        profiles_filtered_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique_fold_filtered.bin',
        compressed_profiles_file='/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/compressed_by_indices_profiles.bin',
        indices_mode=False,
    )

    args = parser.parse_args()

    return args


def test_filtered_profiles(args):
    original_profiles_array = IO.unpack_profiles_file(args.profiles_full_file,
                                                      args.indices_mode,
                                                      do_print = True)
    with open(args.profiles_filtered_file, 'rb') as rf:
        bitstring = rf.read()
    filtered_profiles_array = IO.decompress_profiles_indices(bitstring)

    print(original_profiles_array.shape)
    print(original_profiles_array)
    print(original_profiles_array[6:16,].sum())
    print(filtered_profiles_array.shape)
    print(filtered_profiles_array)
    print(filtered_profiles_array.sum())




def main():
    args = handler()
    #test_bins_fasta(args)
    test_filtered_profiles(args)


if __name__ == "__main__":
    main()


