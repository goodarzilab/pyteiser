import numpy as np
import argparse

import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import IO
import matchmaker
import type_conversions


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--profiles_bin_file", type=str)
    parser.add_argument("--exp_mask_file", help="", type=str)

    parser.set_defaults(
        profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/test_motifs_101.bin",
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin'
    )

    args = parser.parse_args()

    return args


def unpack_profiles_and_mask(args, do_print=False):
    with open(args.profiles_bin_file, 'rb') as rf:
        bitstring = rf.read()
        decompressed_profiles_array = IO.decompress_profiles(bitstring)
        if do_print:
            print("%d profiles have been loaded" % len(decompressed_profiles_array))

    with open(args.exp_mask_file, 'rb') as rf:
        bitstring = rf.read()
        index_array, values_array = IO.decompress_exp_mask_file(bitstring)
        if do_print:
            print("Expression values are provided for %d out of %d transcripts in the reference transcriptome" %
                  (index_array.sum(), index_array.shape[0]))

    try:
        assert (decompressed_profiles_array[0].shape[0] == index_array.shape[0])
    except AssertionError:
        print("Error: occurence profiles were calculated for some other reference transcriptome. The length of the "
              "profiles is %d and the length of the transcriptome provided is %d" %
              (decompressed_profiles_array[0].shape[0], index_array.shape[0]))
        sys.exit(1)

    return decompressed_profiles_array, index_array, values_array



def main():
    args = handler()
    decompressed_profiles_array, index_array, values_array = unpack_profiles_and_mask(args, do_print=False)


if __name__ == "__main__":
    main()
