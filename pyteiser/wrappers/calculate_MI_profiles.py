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
import MI
import matchmaker
import type_conversions


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--profiles_bin_file", type=str)
    parser.add_argument("--exp_mask_file", help="", type=str)
    parser.add_argument("--nbins", help="", type=int)

    parser.set_defaults(
        profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/test_motifs_101.bin",
        #profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/profiles_4-7_4-9_4-6_14-20_30k_1.bin",
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',
        nbins = 3,
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


def quantize_values(index_array, values_array, nbins):
    active_values_array = values_array[index_array]
    quant_values_array = MI.discretize(active_values_array, bins=nbins, disc = "equalfreq")
    return quant_values_array



def calculate_MI_wrapper(decompressed_profiles_array, index_array, quant_values_array):
    for profile in decompressed_profiles_array:
        active_profile = profile[index_array]
        current_MI = MI.mut_info(active_profile, quant_values_array)
        print(current_MI)





def main():
    args = handler()
    decompressed_profiles_array, index_array, values_array = unpack_profiles_and_mask(args, do_print=True)
    quant_values_array = quantize_values(index_array, values_array, args.nbins)
    calculate_MI_wrapper(decompressed_profiles_array, index_array, quant_values_array)

    # proceed with line 179 in mi_find_seed.c



if __name__ == "__main__":
    main()
