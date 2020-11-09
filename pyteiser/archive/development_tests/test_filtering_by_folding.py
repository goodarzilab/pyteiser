import argparse
import numpy as np
import timeit
import hashlib
import struct
import numba

import os
import sys

# to make sure relative import works in order to import test data
current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)
os.chdir(package_home_path)

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
        rna_fastafile='tests/data/test_seqs.fa',
        rna_bin_file='tests/data/test_seqs.bin',
        profiles_full_file='tests/data/profiles.bin',
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


