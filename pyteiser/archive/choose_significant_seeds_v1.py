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

# This script implements greedy search of threshold for statistically significant seeds
# "Greedy" in this context means that it's using a simple objective function
# so "greedy" refers to computational time. Therefore, it's not robust and can find quite
# different (and therefore imprecise) thresholds on each run. The advantage is that it works faster
# For precise threshold identification, something more advanced (like maybe simulated annealing)
# has to be implemented. However, Hani claims that in practice it doesn't matter because seeds are
# redundant and if a good seed didn't pass for some reason there is always a similar seed that did


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--profiles_bin_file", help="file with occurence profiles", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)
    parser.add_argument("--threshold_file", help="file where the threshold ", type=str)

    parser.add_argument("--n_permutations", help="number of permutations for the rnak test for a seed", type=int)
    parser.add_argument("--max_pvalue", help="maximal acceptable p-value", type=int)
    parser.add_argument("--min_zscore", help="maximal acceptable p-value", type=int)
    parser.add_argument("--fastthreshold_jump", help="how many seeds to move down the list in the fast search stage", type=int)
    parser.add_argument("--min_interval", help="when to stop searching in the decreasing intervals phase", type=int)
    parser.add_argument("--min_consecutive_not_passed", help="how many seeds should not pass consecutively "
                                                             "for the search to stop", type=int)

    parser.set_defaults(
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',

        #profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/test_motifs_101.bin",
        profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/profiles_4-7_4-9_4-6_14-20_30k_1.bin",

        #MI_values_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_values/MI_test_motifs_101.bin',
        MI_values_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_values/MI_profiles_4-7_4-9_4-6_14-20_30k_1.bin',

        threshold_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_significancy_threshold/MI_profiles_4-7_4-9_4-6_14-20_30k_1_threshold.bin',

        n_permutations = 1000, # takes 1 second per 100 permutations, Hani's default number of permutations is 1*10^6
        max_pvalue = 0.001, # Hani's default threshold is 1*10^-7
        min_zscore = -1,
        fastthreshold_jump = 50, # Hani's default threshold is 200
        min_interval = 10,
        min_consecutive_not_passed = 10,
    )

    args = parser.parse_args()

    return args


def get_current_statistics(index, MI_values_array, profiles_array,
                           index_array, discr_exp_profile, args):
    profile = profiles_array[index]
    active_profile = profile[index_array]
    current_MI = MI_values_array[index]

    if current_MI == -1:
        return args.max_pvalue + 0.1, args.min_zscore - 0.1

    assert (np.isclose(current_MI, MI.mut_info(active_profile, discr_exp_profile), rtol=1e-10))

    pvalue, z_score = statistic_tests.MI_get_pvalue_and_zscore(active_profile, discr_exp_profile,
                                                               current_MI, args.n_permutations)
    return pvalue, z_score


def determine_thresh_lower_limit(MI_values_array, seed_indices_sorted, seed_pass,
                                 discr_exp_profile, profiles_array, index_array,
                                args, do_print = False):
    last_positive_seed = -1

    counter = 0
    everyone_passed = True

    while (counter < len(seed_indices_sorted)) and everyone_passed:
        index = seed_indices_sorted[counter]
        pvalue, z_score = get_current_statistics(index, MI_values_array, profiles_array,
                                                 index_array, discr_exp_profile, args)

        if pvalue <= args.max_pvalue and z_score >= args.min_zscore:
            # seed passed
            last_positive_seed = counter
            seed_pass[index] = 1
            if do_print:
                print("Seed number %d passed (p=%.5f, z=%.2f)" % (counter, pvalue, z_score))
        else:
            seed_pass[index] = -1
            if do_print:
                print("Seed number %d didn't pass (p=%.5f, z=%.2f)" % (counter, pvalue, z_score))
            everyone_passed = False

        counter += args.fastthreshold_jump

    return last_positive_seed, seed_pass


def decreasing_intervals(last_positive_seed, MI_values_array, seed_indices_sorted,
                         profiles_array, index_array, discr_exp_profile, seed_pass,
                         do_print, args):

    if last_positive_seed >= 0:
        upper_boundary = last_positive_seed # upper limit: last positive seed
        lower_boundary = min(len(MI_values_array) - 1,
                          last_positive_seed + args.fastthreshold_jump) # lower limit
                                # the first negative seed - which is last positive seed + fastthreshold_jump
                                # or, if there were none, the end of the profiles list

        while (lower_boundary - upper_boundary) > args.min_interval:
            counter = upper_boundary + (lower_boundary - upper_boundary) // 2
            index = seed_indices_sorted[counter]
            pvalue, z_score = get_current_statistics(index, MI_values_array, profiles_array,
                                                     index_array, discr_exp_profile, args)
            if pvalue <= args.max_pvalue and z_score >= args.min_zscore:
                # seed passed, go down half interval
                upper_boundary = counter
                last_positive_seed = counter
                seed_pass[index] = 1 # write down that it passed
                if do_print:
                    print("Seed number %d passed (p=%.5f, z=%.2f)" % (counter, pvalue, z_score))
            else:
                # seed didn't pass, go up half interval
                lower_boundary = counter
                seed_pass[index] = -1 # write down that it didn't pass
                if do_print:
                    print("Seed number %d didn't pass (p=%.5f, z=%.2f)" % (counter, pvalue, z_score))
    return last_positive_seed, seed_pass



def search_consec_not_passing_seeds(last_positive_seed, MI_values_array, seed_indices_sorted,
                         profiles_array, index_array, discr_exp_profile, seed_pass,
                         do_print, args):
    if last_positive_seed >= 0:
        number_previous_bad = 0
    else:
        number_previous_bad = 1

    counter = min(len(MI_values_array) - 1, last_positive_seed + args.min_consecutive_not_passed) # why????

    while counter < len(MI_values_array) and number_previous_bad < args.min_consecutive_not_passed:
        index = seed_indices_sorted[counter]

        if seed_pass[index] != 0: # if we have checked this seed before
            check = seed_pass[index]
        else:
            pvalue, z_score = get_current_statistics(index, MI_values_array, profiles_array,
                                                     index_array, discr_exp_profile, args)
            if pvalue <= args.max_pvalue and z_score >= args.min_zscore:
                check = 1 # seed passed
            else:
                check = -1 # seed didn't pass

        if check == 1: # if seed passed
            last_positive_seed = counter # write down last one that passed
            number_previous_bad = 0 # reset number of seeds that didn't pass
            seed_pass[index] = 1 # store info about this seed
            counter += args.min_consecutive_not_passed + 1 # jump several seeds down (because the current one is good)
            # add 1 because 1 is going to be removed immediately below
            if do_print:
                print("Seed number %d passed" % (counter))

        else: # if seed didn't pass
            number_previous_bad += 1
            seed_pass[index] = -1  # store info about this seed
            if do_print:
                print("Seed number %d didn't pass" % (counter))

        counter -= 1

    return last_positive_seed, seed_pass


def determine_mi_threshold(MI_values_array, discr_exp_profile,
                           profiles_array, index_array,
                           args, do_print = False):

    seed_indices_sorted = np.argsort(MI_values_array)[::-1]

    seed_pass = np.zeros(MI_values_array.shape[0], dtype=np.bool) # zero means no info

    if do_print:
        print("Find the lower boundary for the threshold")
    last_positive_seed, seed_pass = determine_thresh_lower_limit(MI_values_array, seed_indices_sorted, seed_pass,
                                                     discr_exp_profile, profiles_array, index_array,
                                                    args, do_print)
    if do_print:
        print("The last seed that passed is: ", last_positive_seed, '\n')
        print("Decreasing intervals phase")
    last_positive_seed, seed_pass = decreasing_intervals(last_positive_seed, MI_values_array, seed_indices_sorted,
                                                     profiles_array, index_array, discr_exp_profile,
                                                     seed_pass, do_print, args)
    if do_print:
        print("The last seed that passed is: ", last_positive_seed, '\n')
        print("Find 10 consecutive seeds that don't pass")
    last_positive_seed, seed_pass = search_consec_not_passing_seeds(last_positive_seed, MI_values_array,
                                                    seed_indices_sorted, profiles_array, index_array,
                                                    discr_exp_profile, seed_pass, do_print, args)
    if do_print:
        print("The last seed that passed is: ", last_positive_seed, '\n')

    IO.write_seed_significancy_threshold(last_positive_seed, args.threshold_file)






def main():
    args = handler()

    # read occurence profiles and expression profile
    profiles_array, index_array, values_array = IO.unpack_profiles_and_mask(args, do_print=False)

    # read precalculated MI values
    MI_values_array, nbins = IO.read_MI_values(args.MI_values_file)

    # find the threshold
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins)
    determine_mi_threshold(MI_values_array, discr_exp_profile,
                           profiles_array, index_array,
                           args, do_print = True)


if __name__ == "__main__":
    main()