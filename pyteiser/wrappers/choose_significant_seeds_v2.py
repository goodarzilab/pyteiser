import numpy as np
import argparse
import numba
import timeit
import math

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

# The original program works in this way
# First, it keeps jumping down taking big steps (fastthreshold_jump) until it finds the first seed that doesn't pass
# Second, it searches for the first seed not to pass in the last jump by decreasing intervals
# Third, it goes down until it finds 10 consecutive seeds that didn't pass
# There is a big fundamental problem with this design: the steps 1 and 2 look for a threshold that defines where the
# first not-passing seed is (point 1). However, the third step looks for a threshold that defines where the last passing
# seed is (or rather where all the seeds are non-passing, which is the same thing) (point 2), and this is a very different
# spot that is located way down the list comparing to the spot where the first non-passing seed is.
# The step 3 checks all the seeds consecutively between point 1 and point 2, which takes super long time and neglects
# the effects of steps 1 and 2 being fast
# Here, in the version 2, this fundamental problem will be fixed: all the 3 steps will be looking for point 2.
# At the steps 1 and 2 I will (1) do bigger jumps and (2) instead of checking 1 seed after each jump, I'll check
# several consecutive seeds (let's say 5 seeds) after each jump, and we stop only if a certain fraction of non-
# passing seeds (say, 4/5) was reached. Then, at stage 2, we do decreasing interval search, also checking several
# seeds at once and looking at the ratio of the ones that passed. At stage 3, we increase the required fraction by a tiny
# bit and we go down until such fraction has been reached
# additionally, since Hani's minimal pvalue is set up in the way that it doesn't pass the threshold even if only 1
# shuffled array gets a pvalue bigger than the one listed, I can (1) remove the threshold parameter from the script and
# (2) stop shuffling as soon as a single MI of bigger value than listed gets encountered. Such hack will speed up the
# threshold search


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--profiles_bin_file", help="file with occurence profiles", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)
    parser.add_argument("--threshold_file", help="file where the threshold ", type=str)

    parser.add_argument("--n_permutations", help="number of permutations for the rnak test for a seed", type=int)
    parser.add_argument("--max_pvalue", help="maximal acceptable p-value", type=int)
    parser.add_argument("--min_zscore", help="maximal acceptable p-value", type=int)
    parser.add_argument("--step_1_jump", help="", type=int)
    parser.add_argument("--step_2_min_interval", help="", type=int)
    parser.add_argument("--step_1_min_fraction", help="", type=float)
    parser.add_argument("--step_2_min_fraction", help="", type=float)
    parser.add_argument("--step_3_min_fraction", help="", type=float)

    parser.set_defaults(
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',

        #profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/test_motifs_101.bin",
        profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/profiles_4-7_4-9_4-6_14-20_30k_1.bin",

        #MI_values_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_values/MI_test_motifs_101.bin',
        MI_values_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_values/MI_profiles_4-7_4-9_4-6_14-20_30k_1.bin',

        threshold_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_significancy_threshold/MI_profiles_4-7_4-9_4-6_14-20_30k_1_threshold.bin',

        n_permutations = 1000, # takes 1 second per 100 permutations, Hani's default number of permutations is 1*10^6
        max_pvalue = 0.0001, # Hani's default threshold is 1*10^-7
        min_zscore = -1,

        step_1_jump = 100, # Hani's default jump is 200
        step_2_min_interval=10,

        step_1_min_fraction = 0.8,
        step_2_min_fraction= 0.8,
        step_3_min_fraction= 0.9,
    )

    args = parser.parse_args()

    return args

# let's say I want
# instead of defining two separate parameters: (1) how many

def find_fraction_denominator(x):
    rounded_denominator = round(1 / (1 - x), 4)
    return math.ceil(rounded_denominator)



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


def check_one_seed(index, counter, seed_pass, MI_values_array,
                   profiles_array, index_array, discr_exp_profile,
                   args, do_print):
    pv_defined = False
    if seed_pass[index] != 0:  # if we have checked this seed before
        check = seed_pass[index]
    else:
        pvalue, z_score = get_current_statistics(index, MI_values_array, profiles_array,
                                                 index_array, discr_exp_profile, args)
        pv_defined = True
        if pvalue <= args.max_pvalue and z_score >= args.min_zscore:
            check = 1  # seed passed
        else:
            check = -1  # seed didn't pass

    if check == 1:
        did_pass = True
        seed_pass[index] = 1
        if do_print:
            if pv_defined:
                print("Seed number %d passed (p=%.5f, z=%.2f)" % (counter, pvalue, z_score))
            else:
                print("Seed number %d passed" % (counter))
    elif check == -1:
        did_pass = False
        seed_pass[index] = -1
        if do_print:
            if pv_defined:
                print("Seed number %d didn't pass (p=%.5f, z=%.2f)" % (counter, pvalue, z_score))
            else:
                print("Seed number %d didn't pass" % (counter))
    else:
        print("Error!")
        sys.exit(1)

    return did_pass, seed_pass


def check_N_consecutive(minimal_fraction, start_point,
        MI_values_array, seed_indices_sorted, seed_pass,
                                 discr_exp_profile, profiles_array, index_array,
                                args, do_print = False):
    number_not_passed = 0
    denominator = find_fraction_denominator(minimal_fraction)
    go_until =  start_point + denominator
    if go_until > len(MI_values_array) - 1:
        print("There aren't enough seeds to check")
        return False, seed_pass

    for counter in range(start_point, start_point + denominator):
        index = seed_indices_sorted[counter]
        did_seed_pass, seed_pass = check_one_seed(index, counter, seed_pass, MI_values_array,
                                       profiles_array, index_array, discr_exp_profile,
                                       args, do_print)
        if not did_seed_pass:
            # seed didn't pass, write it down to number_not_passed variable
            number_not_passed += 1

    fraction_not_passed = float(number_not_passed) / denominator
    if fraction_not_passed >= minimal_fraction:
        is_over_threshold = True
    else:
        is_over_threshold = False

    return is_over_threshold, seed_pass


def step_1_determine_thresh_lower_limit(MI_values_array, seed_indices_sorted, seed_pass,
                                 discr_exp_profile, profiles_array, index_array,
                                args, do_print = False):
    last_positive_seed = -1

    counter = 0
    is_over_threshold = False

    while (counter < len(seed_indices_sorted)) and not is_over_threshold:
        is_over_threshold, seed_pass = check_N_consecutive(args.step_1_min_fraction,
                            counter, MI_values_array, seed_indices_sorted, seed_pass,
                            discr_exp_profile, profiles_array, index_array,
                            args, do_print)
        last_positive_seed = counter
        counter += args.step_1_jump

    if not is_over_threshold:
        print("No lower boundary for threshold found!")
        return len(MI_values_array) - 1, seed_pass

    return last_positive_seed, seed_pass


def step_2_decreasing_intervals(last_positive_seed, MI_values_array, seed_indices_sorted,
                         profiles_array, index_array, discr_exp_profile, seed_pass,
                         do_print, args):

    if last_positive_seed >= 0:
        upper_boundary = last_positive_seed - args.step_1_jump # upper limit: last positive seed - one jump (just in case)
        lower_boundary = min(len(MI_values_array) - 1,
                          last_positive_seed + args.step_1_jump) # lower limit
                                # the first negative seed - which is last positive seed + step_1_jump
                                # or, if there were none, the end of the profiles list

        while (lower_boundary - upper_boundary) > args.step_2_min_interval:
            counter = upper_boundary + (lower_boundary - upper_boundary) // 2
            is_over_threshold, seed_pass = check_N_consecutive(args.step_1_min_fraction,
                                                               counter, MI_values_array, seed_indices_sorted, seed_pass,
                                                               discr_exp_profile, profiles_array, index_array,
                                                               args, do_print)
            if not is_over_threshold:
                # some seed passed, go down half interval
                upper_boundary = counter
                last_positive_seed = counter
            else:
                # enough seeds didn't pass, go up half interval
                lower_boundary = counter

    return last_positive_seed, seed_pass


def get_current_fraction_from_array(x):
    number_not_passed = np.logical_not(x).sum()
    return float(number_not_passed) / x.shape[0]


def step_3_confirm_consec_not_passing_seeds(last_positive_seed, MI_values_array, seed_indices_sorted,
                         profiles_array, index_array, discr_exp_profile, seed_pass,
                         do_print, args):
    if last_positive_seed < 0: # if even the first seed didn't pass in the 1st step, the last_positive_seed variable
        last_positive_seed = 0 # holds -1, but we need it to hold 0

    denominator = find_fraction_denominator(args.step_3_min_fraction)
    if last_positive_seed + denominator > len(MI_values_array):
        print("Reached the end of file, no threshold found! Suppose that all the seeds are significant")
        sys.exit(1)

    current_passed_array = np.zeros(denominator, dtype=np.bool)
    for i in range(denominator):
        counter = i + last_positive_seed + 1
        index = seed_indices_sorted[counter]
        did_seed_pass, seed_pass = check_one_seed(index, counter, seed_pass, MI_values_array,
                                       profiles_array, index_array, discr_exp_profile,
                                       args, do_print)
        current_passed_array[i] = did_seed_pass

    current_fraction = get_current_fraction_from_array(current_passed_array)
    fraction_exceeded = (current_fraction > args.step_3_min_fraction)

    counter = last_positive_seed + denominator

    while counter < len(MI_values_array) and not fraction_exceeded:
        index = seed_indices_sorted[counter]
        did_seed_pass, seed_pass = check_one_seed(index, counter, seed_pass, MI_values_array,
                                                  profiles_array, index_array, discr_exp_profile,
                                                  args, do_print)
        current_passed_array = np.roll(current_passed_array, -1)
        current_passed_array[-1] = did_seed_pass
        current_fraction = get_current_fraction_from_array(current_passed_array)
        fraction_exceeded = (current_fraction > args.step_3_min_fraction)
        counter += 1

    if not fraction_exceeded:
        print("Reached the end of file, no threshold found! Suppose that all the seeds are significant")
        sys.exit(1)

    last_positive_seed = counter - denominator - 1
    return last_positive_seed, seed_pass


def determine_mi_threshold(MI_values_array, discr_exp_profile,
                           profiles_array, index_array,
                           args, do_print = False):

    seed_indices_sorted = np.argsort(MI_values_array)[::-1]

    seed_pass = np.zeros(MI_values_array.shape[0], dtype=np.int8) # zero means no info

    if do_print:
        print("Find the lower boundary for the threshold")
    last_positive_seed, seed_pass = step_1_determine_thresh_lower_limit(MI_values_array, seed_indices_sorted, seed_pass,
                                                     discr_exp_profile, profiles_array, index_array,
                                                    args, do_print)
    if do_print:
        print("The last seed that passed is: ", last_positive_seed, '\n')
        print("Decreasing intervals phase")
    last_positive_seed, seed_pass = step_2_decreasing_intervals(last_positive_seed, MI_values_array, seed_indices_sorted,
                                                     profiles_array, index_array, discr_exp_profile,
                                                     seed_pass, do_print, args)
    if do_print:
        print("The last seed that passed is: ", last_positive_seed, '\n')
        print("Find 10 consecutive seeds that don't pass")

    last_positive_seed, seed_pass = step_3_confirm_consec_not_passing_seeds(last_positive_seed, MI_values_array,
                                            seed_indices_sorted, profiles_array, index_array, discr_exp_profile,
                                            seed_pass, do_print, args)

    if do_print:
        print("The last seed that passed is: ", last_positive_seed, '\n')

    IO.write_seed_significancy_threshold(last_positive_seed, args.threshold_file)


def main():
    args = handler()

    # read occurence profiles and expression profile
    profiles_array, index_array, values_array = IO.unpack_profiles_and_mask(args.profiles_bin_file,
                                                                            args.exp_mask_file, do_print=False)
    # read precalculated MI values
    MI_values_array, nbins = IO.read_MI_values(args.MI_values_file)

    # find the threshold
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins)
    determine_mi_threshold(MI_values_array, discr_exp_profile,
                           profiles_array, index_array,
                           args, do_print = True)


if __name__ == "__main__":
    main()
