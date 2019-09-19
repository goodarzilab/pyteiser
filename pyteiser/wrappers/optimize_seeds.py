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
    parser.add_argument("--combined_passed_seeds_filename", help="file with the seeds that have passed the threshold", type=str)
    parser.add_argument("--rna_bin_file", help="referense transcriptome in binary format", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)

    parser.add_argument("--maxfreq", help="", type=float)
    parser.add_argument("--n_permutations", help="number of permutations for the rank test for a seed", type=int)
    parser.add_argument("--jackknife_n_samples", help="how many permutations to do in jackknife test", type=int)
    parser.add_argument("--jackknife_fraction_retain", help="what fraction of the sample to retain for each test",
                                                                                    type=float)
    parser.add_argument("--jackknife_min_fraction_passed", help="what fraction of all iterations should"
                                                                "pass to consider the motif robust", type=float)

    parser.set_defaults(
        combined_passed_seeds_filename='/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/seeds_passed_100k_tarbp2_utrs.bin',
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',

        maxfreq = 0.5, # default value from Hani's program is 0.5
        n_permutations=1000,  # takes 1 second per 100 permutations, Hani's default number of permutations is 1*10^6
        jackknife_n_samples = 10,
        jackknife_fraction_retain = 0.66,
        jackknife_min_fraction_passed = 0.6,

    )

    args = parser.parse_args()

    return args



def are_there_better_motifs(n_modified_motifs, n_seqs_list, discr_exp_profile, nbins,
                            bestmi, n_bestmotif, lastmyfreq, args):

    for j in range(len(n_modified_motifs)):
        # limit it to the expressed sequences only
        # change the matchmaking procedure to incorporate degenerative nucleotides
        current_profile, time_spent = matchmaker.calculate_profile_one_motif(n_modified_motifs[j],
                                                                             n_seqs_list,
                                                                            is_degenerate = True)
        myfreq = current_profile.sum() / float(len(n_seqs_list))
        tempmi = MI.mut_info(current_profile, discr_exp_profile, x_bins=2, y_bins=nbins)

        if tempmi > bestmi and current_profile.sum() > 10 and (myfreq < args.maxfreq or myfreq < lastmyfreq):
            n_bestmotif = n_modified_motifs[j].copy()
            w_bestmotif = type_conversions.n_to_w_sequence(n_bestmotif)
            bestmi = tempmi
            lastmyfreq = myfreq
            print("new motif (mi = %.4f): %s\n", bestmi, w_bestmotif.print_sequence(return_string=True))
    return bestmi, lastmyfreq, n_bestmotif


# optimize sequence of all the positions individually in random order
def optimize_motif_sequence(counter, index, n_motifs_list, n_seqs_list,
                            discr_exp_profile, nbins, active_profile,
                            init_best_mymi, n_bestmotif, lastmyfreq, args):
    # print initial motif
    bestmi = init_best_mymi
    print("Optimzing the sequence of motif %d" % counter)
    w_bestmotif = type_conversions.n_to_w_sequence(n_bestmotif)
    print("initial motif (mi = %.4f): %s", bestmi, w_bestmotif.print_sequence(return_string=True))

    # create a random index so that we optimize each position not from left to right but rather in random order
    k_inc = np.arange(n_motifs_list[index].length)
    k_shu = np.random.permutation(k_inc)

    # optimize motif
    for k in n_motifs_list[index].length:
        position = k_shu[k]
        w_modified_motifs = modify_seed.modify_base(n_motifs_list[index], position)
        n_modified_motifs = type_conversions.w_to_n_motifs_list(w_modified_motifs)
        bestmi, lastmyfreq, n_bestmotif = are_there_better_motifs(n_modified_motifs,
                                                    n_seqs_list, discr_exp_profile, nbins,
                                                    bestmi, n_bestmotif, lastmyfreq, args)
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


def optimize_motifs(number_signigicant_seeds, MI_values_array, discr_exp_profile, nbins,
                    profiles_array, index_array, n_motifs_list, n_seqs_list, args, do_print = False):

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


        n_bestmotif = n_motifs_list[index].copy()

        # initial mi value
        init_best_mymi = MI.mut_info(active_profile, discr_exp_profile, x_bins=2, y_bins=nbins)
        lastmyfreq = active_profile.sum() / float(active_profile.shape[0])

        print("Initial MI = %.5f\n", init_best_mymi)

        bestmi, lastmyfreq, n_bestmotif = optimize_motif_sequence(counter, index,
                                n_motifs_list, n_seqs_list,
                                discr_exp_profile, nbins, active_profile,
                                init_best_mymi, n_bestmotif, lastmyfreq, args)
        bestmi, lastmyfreq, n_bestmotif = elongate_motif(counter, index,
                       n_motifs_list, n_seqs_list,
                       discr_exp_profile, nbins, active_profile,
                       bestmi, n_bestmotif, lastmyfreq, args)

        bestmotif_profile, _ = matchmaker.calculate_profile_one_motif(n_bestmotif, n_seqs_list)
        bestmotif_active_profile = bestmotif_profile[index_array]
        bestmotif_mi = MI.mut_info(bestmotif_active_profile, discr_exp_profile, x_bins=2, y_bins=nbins)
        pvalue, z_score = statistic_tests.MI_get_pvalue_and_zscore(bestmotif_active_profile, discr_exp_profile, nbins,
                                                                   bestmotif_mi, args.n_permutations)
        print("The findal z-score is: %.3f", z_score)

        passed_jacknife = statistic_tests.jackknife_test(args.jackknife_n_samples,
                                                        args.jackknife_fraction_retain,
                                                        args.jackknife_min_fraction_passed)


        pv_defined = True
        if pvalue <= args.max_pvalue and z_score >= args.min_zscore:
            check = 1  # seed passed
        else:
            check = -1  # seed didn't pass




        # TODO: stopped at line 271 of mi_optimize.c









def read_input_files(seeds_filename_full, rna_bin_filename):
    seqs_dict, seqs_order = IO.read_rna_bin_file(rna_bin_filename)
    w_motifs_list = IO.read_motif_file(seeds_filename_full)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]
    n_motifs_list = type_conversions.w_to_n_motifs_list(w_motifs_list)
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)

    return n_motifs_list, n_seqs_list





def main():
    args = handler()

    # read seeds and sequences
    n_motifs_list, n_seqs_list = read_input_files(args.seed_file, args.rna_bin_file)

    # read expression profile
    index_array, values_array = IO.unpack_mask_file(args.exp_mask_file)

    # read precalculated MI values
    MI_values_array, nbins = IO.read_MI_values(args.MI_values_file)

    # read precalculated threshold
    number_signigicant_seeds = IO.read_seed_significancy_threshold(args.threshold_file)


    # optimize motifs
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins)
    optimize_motifs(number_signigicant_seeds,
                    MI_values_array, discr_exp_profile, nbins,
                       profiles_array, index_array,
                    n_motifs_list, n_seqs_list,
                       args, do_print=True)





if __name__ == "__main__":
    main()
