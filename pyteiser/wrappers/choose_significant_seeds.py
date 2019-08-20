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

MASK_OUT_SEED_VALUE = np.float64(-1)


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--profiles_bin_file", help="file with occurence profiles", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)
    parser.add_argument("--MI_values_file", help="file with precalculated MI values for individual profiles", type=str)

    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)
    parser.add_argument("--min_occurences", help="minimal number of seed occurence in the transcriptome"
                                                 " for a seed to be considered", type=int)


    parser.set_defaults(
        profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/test_motifs_101.bin",
        #profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/profiles_4-7_4-9_4-6_14-20_30k_1.bin",
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',
        MI_values_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_values/MI_test_motifs_101.bin',

        nbins = 3,
        min_occurences = 5,
    )

    args = parser.parse_args()

    return args



def determine_mi_threshold(MI_values_array, discr_exp_profile):
    seed_indices_sorted = np.argsort(MI_values_array)[::-1]
    seed_pass = np.zeros(MI_values_array.shape[0], dtype=np.bool)

    for counter, index in enumerate(seed_indices_sorted):
        print(counter, index, MI_values_array[index])
        #seed_max_rank_test(discr_exp_profile)


    # print(MI_values_sorted)

def seed_max_rank_test(discr_exp_profile):
    shuffled_expr = np.random.permutation(discr_exp_profile)


def main():
    args = handler()
    decompressed_profiles_array, index_array, values_array = IO.unpack_profiles_and_mask(args, do_print=False)
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, args.nbins)
    with open(args.MI_values_file, 'rb') as rf:
        bitstring = rf.read()
    MI_values_array = IO.decompres_MI_values(bitstring)

    determine_mi_threshold(MI_values_array, discr_exp_profile)


    # print(MI_values_array.shape)
    # print(MI_values_array[0])
    # print(np.float32(MI_values_array[0]))

    #print(np.array(MI_values_array, dtype=np.float32)*10000)

    # MI_values_array = calculate_MI_for_seeds(decompressed_profiles_array, index_array, discr_exp_profile,
    #                                      args.min_occurences, do_print = True)
    # IO.write_MI_values(MI_values_array, args.MI_values_file)


    # proceed with line 179 in mi_find_seed.c



if __name__ == "__main__":
    main()
