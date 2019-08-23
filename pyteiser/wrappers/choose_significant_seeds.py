import numpy as np
import argparse
import numba


import timeit

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
import statistic_tests

MASK_OUT_SEED_VALUE = np.float64(-1)


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--profiles_bin_file", help="file with occurence profiles", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)

    parser.add_argument("--n_permutations", help="number of permutations for the rnak test for a seed", type=int)
    parser.add_argument("--max_pvalue", help="maximal acceptable p-value", type=int)
    parser.add_argument("--min_zscore", help="maximal acceptable p-value", type=int)
    parser.add_argument("--fastthreshold_jump", help="how many seeds to move down the list in the fast search stage", type=int)

    parser.set_defaults(
        profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/test_motifs_101.bin",
        #profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/profiles_4-7_4-9_4-6_14-20_30k_1.bin",
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',
        MI_values_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_values/MI_test_motifs_101.bin',

        n_permutations = 100, # takes 1 second per 100 permutations
        max_pvalue = 0.01, # Hani's default threshold is 0.0000001
        min_zscore = -1,
        fastthreshold_jump = 1, # Hani's default threshold is 200
    )

    args = parser.parse_args()

    return args



def determine_mi_threshold(MI_values_array, discr_exp_profile,
                           profiles_array, index_array,
                           args, do_print = False):
    seed_indices_sorted = np.argsort(MI_values_array)[::-1]
    # seed_pass contains 1-pvalue values for seeds as calculated by a permutation test in seed_max_rank_test
    seed_pass = np.zeros(MI_values_array.shape[0], dtype=np.bool) # zero means no info
    last_positive_seed = -1

    for counter in range(0, len(seed_indices_sorted), args.fastthreshold_jump):
        index = seed_indices_sorted[counter]
    #for counter, index in enumerate(seed_indices_sorted):
        profile = profiles_array[index]
        active_profile = profile[index_array]
        current_MI = MI_values_array[index]

        if current_MI == -1:
            seed_pass[index] = 0
            continue

        assert(np.isclose(current_MI, MI.mut_info(active_profile, discr_exp_profile), rtol=1e-10))

        pvalue, z_score = statistic_tests.MI_get_pvalue_and_zscore(active_profile, discr_exp_profile,
                           current_MI, args.n_permutations)
        # time it
        # time_norm = timeit.timeit(lambda: statistic_tests.MI_get_pvalue_and_zscore(active_profile, discr_exp_profile,
        #                    current_MI, n_permutations), number=2)
        # print(time_norm)

        if pvalue > args.max_pvalue or z_score < args.min_zscore:
            seed_pass[index] = 1
            if do_print:
                print("Seed number %d didn't pass (p=%.5f, z=%.2f" % (counter, pvalue, z_score))
            break
        else:
            last_positive_seed = counter
            seed_pass[index] = -1
            if do_print:
                print("Seed number %d passed (p=%.5f, z=%.2f" % (counter, pvalue, z_score))

        # if counter > 2:
        #     break




def main():
    args = handler()

    # read occurence profiles and expression profile
    profiles_array, index_array, values_array = IO.unpack_profiles_and_mask(args, do_print=False)

    # read precalculated MI values
    with open(args.MI_values_file, 'rb') as rf:
        bitstring = rf.read()
    MI_values_array, nbins = IO.decompres_MI_values(bitstring)

    # find the threshold
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins)
    determine_mi_threshold(MI_values_array, discr_exp_profile,
                           profiles_array, index_array,
                           args, do_print = True)


    # proceed with line 179 in mi_find_seed.c



if __name__ == "__main__":
    main()
