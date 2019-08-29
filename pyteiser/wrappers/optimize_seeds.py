import numpy as np
import argparse
import hashlib
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
import structures
import statistic_tests

MASK_OUT_SEED_VALUE = np.float64(-1)

NUMBER_MODIFIED_MOTIFS_1 = 15
NUMBER_MODIFIED_MOTIFS_2 = 46


# This script implements greedy search of threshold for statistically significant seeds
# "Greedy" in this context means that it's using a simple objective function
# so "greedy" refers to computational time. Therefore, it's not robust and can find quite
# different (and therefore imprecise) thresholds on each run. The advantage is that it works faster
# For precise threshold identification, something more advanced (like maybe simulated annealing)
# has to be implemented. However, Hani claims that in practice it doesn't matter because seeds are
# redundant and if a good seed didn't pass for some reason there is always a similar seed that did



# add seeds file to the parameters!


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

        n_permutations = 100, # takes 1 second per 100 permutations
        max_pvalue = 0.01, # Hani's default threshold is 0.0000001
        min_zscore = -1,
        fastthreshold_jump = 50, # Hani's default threshold is 200
        min_interval = 10,
        min_consecutive_not_passed = 10,
    )

    args = parser.parse_args()

    return args



def optimize_motifs(number_signigicant_seeds, MI_values_array, discr_exp_profile,
                    profiles_array, index_array, args, do_print = False):

    seed_indices_sorted = np.argsort(MI_values_array)[::-1]
    signif_indices = seed_indices_sorted[0 : number_signigicant_seeds]

    print("There are %d signif seeds in total" % number_signigicant_seeds)
    for counter, index in enumerate(signif_indices):
        profile = profiles_array[index]
        active_profile = profile[index_array]
        current_MI = MI_values_array[index]

        # check how much information it adds to the previous guys
        if opt_count > 0 and minr>0:
            minratio = minCondInfoNormalized()
        else:
            minratio = minr + 1

        if minratio < minr:
            print("seed %d killed by motif %d (ratio=%f).\n", i, midx, minratio)
            continue
        else:
            print("optimizing")

        if doonlypositive:
            r = pearson_int(M_q, E_q, seq_count)
            if r < 0:
                print("seed %d killed due to negative association (pearson=%4.3f)\n", i, r)
                continue

        lastmyfreq = hits / seq_count
        best_lastmyfreq = lastmyfreq

        if do_optimize:
            # initial mi value
            init_best_mymi = MI.mut_info(active_profile, discr_exp_profile)
            print("Initial MI = %.5f\n", init_best_mymi)

            # create a random index
            k_inc = np.arange(motifs[index].length)
            k_shu = np.random.permutation(k_inc)

            # optimize motif
            bestmi = init_best_mymi
            print("Optimzing the sequence of motif %d" % counter)
            print("initial motif (mi = %.4f): %s\n", bestmi, bestmotif.print_sequence(do_return=True))

            for k in motifs[index].length:
                pos = k_shu[k]
                modified_motifs = [0] * NUMBER_MODIFIED_MOTIFS_1
                for j in range(NUMBER_MODIFIED_MOTIFS_1):
                    modified_motifs[j] = structures.w_motif(? ?)
                    





def main():
    args = handler()

    #read occurence profiles and expression profile
    profiles_array, index_array, values_array = IO.unpack_profiles_and_mask(args, do_print=False)

    # read precalculated MI values
    MI_values_array, nbins = IO.read_MI_values(args.MI_values_file)

    # read precalculated threshold
    number_signigicant_seeds = IO.read_seed_significancy_threshold(args.threshold_file)


    # optimize motifs
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins)
    optimize_motifs(number_signigicant_seeds,
                    MI_values_array, discr_exp_profile,
                       profiles_array, index_array,
                       args, do_print=True)





if __name__ == "__main__":
    main()
