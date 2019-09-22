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
import glob_var
import modify_seed
import type_conversions
import matchmaker
import statistic_tests

MASK_OUT_SEED_VALUE = np.float64(-1)

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

    parser.add_argument("--unique_seeds_filename", help="output: best representatives of each family", type=str)
    parser.add_argument("--unique_profiles_filename", help="output: profiles of best representatives of each family",
                                                                                                        type=str)
    parser.add_argument("--families_classification_filename", help="output: classification of all the passed seeds"
                                                                   "to unique families", type=str)




    parser.add_argument("--rna_bin_file", help="referense transcriptome in binary format", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)

    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)
    parser.add_argument("--maxfreq", help="", type=float)
    parser.add_argument("--n_permutations", help="number of permutations for the rank test for a seed", type=int)
    parser.add_argument("--jackknife_n_samples", help="how many permutations to do in jackknife test", type=int)
    parser.add_argument("--jackknife_fraction_retain", help="what fraction of the sample to retain for each test",
                                                                                    type=float)
    parser.add_argument("--jackknife_min_fraction_passed", help="what fraction of all iterations should"
                                                                "pass to consider the motif robust", type=float)

    parser.set_defaults(
        # unique_seeds_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/seeds_unique_100k_tarbp2_utrs.bin",
        # unique_profiles_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/profiles_unique_100k_tarbp2_utrs.bin",
        # families_classification_filename="/Users/student/Documents/hani/programs/pyteiser/data/seeds_family_classification/seeds_4-7_4-9_4-6_14-20_combined/seeds_unique_100k_tarbp2_utrs_classification.bin",


        # unique_seeds_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/seeds_unique_100k_snrnpa1.bin",
        # unique_profiles_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/profiles_unique_100k_snrnpa1.bin",
        # families_classification_filename="/Users/student/Documents/hani/programs/pyteiser/data/seeds_family_classification/seeds_4-7_4-9_4-6_14-20_combined/seeds_unique_100k_snrnpa1_classification.bin",


        unique_seeds_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds_unique.bin",
        unique_profiles_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique.bin",
        families_classification_filename="/Users/student/Documents/hani/programs/pyteiser/data/seeds_family_classification/seeds_4-7_4-9_4-6_14-20_combined/test_1_2_classification.bin",

        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',

        nbins=15,
        maxfreq = 0.5, # default value from Hani's program is 0.5
        n_permutations=1000,  # takes 1 second per 100 permutations, Hani's default number of permutations is 1*10^6
        jackknife_n_samples = 10,
        jackknife_fraction_retain = 0.66,
        jackknife_min_fraction_passed = 0.6,

    )

    args = parser.parse_args()

    return args


def import_modules():
    package_home_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if package_home_path not in sys.path:
        sys.path.append(package_home_path)

    global MI
    global IO

    import MI
    import IO



def are_there_better_motifs(n_modified_motifs, seqs_of_interest, discr_exp_profile, nbins,
                            bestmi, n_bestmotif, lastmyfreq, args, do_print = True):

    for j in range(len(n_modified_motifs)):
        # limit it to the expressed sequences only
        # change the matchmaking procedure to incorporate degenerative nucleotides
        current_profile, time_spent = matchmaker.calculate_profile_one_motif(n_modified_motifs[j],
                                                                             seqs_of_interest,
                                                                            is_degenerate = True)
        myfreq = current_profile.values.sum() / float(len(seqs_of_interest))
        tempmi = MI.mut_info(current_profile.values, discr_exp_profile, x_bins=2, y_bins=nbins)

        print(tempmi)

        if tempmi > bestmi and current_profile.sum() > 10 and (myfreq < args.maxfreq or myfreq < lastmyfreq):
            n_bestmotif = structures.copy_n_motif(n_modified_motifs[j])
            w_bestmotif = type_conversions.n_to_w_motif(n_bestmotif)
            bestmi = tempmi
            lastmyfreq = myfreq
            if do_print:
                print("new motif (mi = %.4f): %s\n" % (bestmi, w_bestmotif.print_sequence(return_string=True)))
    return bestmi, lastmyfreq, n_bestmotif


# optimize sequence of all the positions individually in random order
def optimize_motif_sequence(n_bestmotif, init_best_MI, seqs_of_interest,
                            discr_exp_profile, nbins, lastmyfreq, args,
                            do_print = False):
    bestmi = init_best_MI

    # create a random index so that we optimize each position not from left to right but rather in random order
    k_inc = np.arange(n_bestmotif.length)
    k_shu = np.random.permutation(k_inc)

    # optimize motif
    for k in range(n_bestmotif.length):
        if do_print:
            print("Modifying position ", k+1)
        position = k_shu[k]
        w_modified_motifs = modify_seed.modify_base(n_bestmotif, position)
        n_modified_motifs = type_conversions.w_to_n_motifs_list(w_modified_motifs)
        bestmi, lastmyfreq, n_bestmotif = are_there_better_motifs(n_modified_motifs,
                                                    seqs_of_interest, discr_exp_profile, nbins,
                                                    bestmi, n_bestmotif, lastmyfreq, args,
                                                    do_print = do_print)
    return bestmi, lastmyfreq, n_bestmotif


def elongate_motif(counter, index, n_motifs_list, n_seqs_list,
                    discr_exp_profile, nbins, active_profile,
                    bestmi, n_bestmotif, lastmyfreq, args):
    print("Elongating the motif %d" % counter)

    premi = bestmi

    # TODO: how to I to a first iteration of this?
    # TODO: also why premi >= bestmi, shouldn't it go the other way??

    while premi >= bestmi and n_bestmotif.length > n_motifs_list[index].length:
        n_elongated_motifs = modify_seed.elongate_motif(n_motifs_list[index])
        for j in range(len(n_elongated_motifs)):
            bestmi, lastmyfreq, n_bestmotif = are_there_better_motifs(n_elongated_motifs,
                                                          n_seqs_list, discr_exp_profile, nbins,
                                                          bestmi, n_bestmotif, lastmyfreq, args)
    return bestmi, lastmyfreq, n_bestmotif


# def optimize_motifs(number_signigicant_seeds, discr_exp_profile, nbins,
#                     profiles_array, index_array, n_motifs_list, n_seqs_list, args, do_print = False):
def optimize_motifs(seeds_initial, profiles_initial,
                    discr_exp_profile, nbins, index_array, seqs_of_interest,
                    args, do_print = True):

    print("Starting with %d initial seeds" % len(seeds_initial))
    for i, motif in enumerate(seeds_initial):
        profile = profiles_initial[i]
        active_profile = profile[index_array]

        n_bestmotif = seeds_initial[i].copy()

        # initial mi value
        init_best_MI = MI.mut_info(active_profile, discr_exp_profile, x_bins=2, y_bins=nbins)
        lastmyfreq = active_profile.sum() / float(active_profile.shape[0])

        if do_print:
            print("Optimzing the sequence of motif %d" % i)
            print("Initial MI = %.5f" % init_best_MI)
            w_bestmotif = type_conversions.n_to_w_motif(n_bestmotif)
            print("initial motif: %s" % w_bestmotif.print_sequence(return_string=True))

        bestmi, lastmyfreq, n_bestmotif = optimize_motif_sequence(n_bestmotif, init_best_MI, seqs_of_interest,
                            discr_exp_profile, nbins, lastmyfreq, args, do_print = do_print)


        # bestmi, lastmyfreq, n_bestmotif = elongate_motif(counter, index,
        #                n_motifs_list, n_seqs_list,
        #                discr_exp_profile, nbins, active_profile,
        #                bestmi, n_bestmotif, lastmyfreq, args)
        #
        # bestmotif_profile, _ = matchmaker.calculate_profile_one_motif(n_bestmotif, n_seqs_list)
        # bestmotif_active_profile = bestmotif_profile[index_array]
        # bestmotif_mi = MI.mut_info(bestmotif_active_profile, discr_exp_profile, x_bins=2, y_bins=nbins)
        # pvalue, z_score = statistic_tests.MI_get_pvalue_and_zscore(bestmotif_active_profile, discr_exp_profile, nbins,
        #                                                            bestmotif_mi, args.n_permutations)
        # print("The findal z-score is: %.3f", z_score)
        #
        # passed_jacknife = statistic_tests.jackknife_test(args.jackknife_n_samples,
        #                                                 args.jackknife_fraction_retain,
        #                                                 args.jackknife_min_fraction_passed)
        #
        #
        # pv_defined = True
        # if pvalue <= args.max_pvalue and z_score >= args.min_zscore:
        #     check = 1  # seed passed
        # else:
        #     check = -1  # seed didn't pass

        # TODO: stopped at line 271 of mi_optimize.c

        break









def read_sequences(rna_bin_filename):
    seqs_dict, seqs_order = IO.read_rna_bin_file(rna_bin_filename)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)

    return n_seqs_list





def main():
    import_modules()
    args = handler()

    n_seqs_list = read_sequences(args.rna_bin_file)
    index_array, values_array = IO.unpack_mask_file(args.exp_mask_file)
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins = args.nbins)
    seeds_initial = IO.read_motif_file(args.unique_seeds_filename)
    profiles_initial = IO.unpack_profiles_file(args.unique_profiles_filename)

    seqs_of_interest = [n_seqs_list[x] for x in range(index_array.shape[0]) if index_array[x]]
    optimize_motifs(seeds_initial, profiles_initial,
                    discr_exp_profile, args.nbins, index_array, seqs_of_interest,
                    args, do_print=True)


    # classif_array = IO.read_classification_array(args.families_classification_filename)

    # read sequences
    # n_seqs_list = read_sequences(args.rna_bin_file)

    # read expression profile
    # index_array, values_array = IO.unpack_mask_file(args.exp_mask_file)


    # optimize motifs
    # discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins)
    # optimize_motifs(number_signigicant_seeds,
    #                 MI_values_array, discr_exp_profile, nbins,
    #                    profiles_array, index_array,
    #                 n_motifs_list, n_seqs_list,
    #                    args, do_print=True)





if __name__ == "__main__":
    main()
