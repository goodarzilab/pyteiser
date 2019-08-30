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

    parser.add_argument("--profiles_bin_file", help="file with occurence profiles", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)
    parser.add_argument("--threshold_file", help="file where the threshold ", type=str)


    parser.add_argument("--seed_file", help="file with the seeds corresponding to the profiles", type=str)
    parser.add_argument("--rna_bin_file", help="", type=str)

    parser.add_argument("--maxfreq", help="", type=float)

    parser.set_defaults(
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',

        #profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/test_motifs_101.bin",
        profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/profiles_4-7_4-9_4-6_14-20_30k_1.bin",

        #MI_values_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_values/MI_test_motifs_101.bin',
        MI_values_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_values/MI_profiles_4-7_4-9_4-6_14-20_30k_1.bin',

        threshold_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_significancy_threshold/MI_profiles_4-7_4-9_4-6_14-20_30k_1_threshold.bin',

        seed_file='/Users/student/Documents/hani/programs/pyteiser/data/test_seeds/seeds_4-7_4-9_4-6_14-20_30k_1.bin',
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',

        maxfreq = 0.5, # default value from Hani's program is 0.5
    )

    args = parser.parse_args()

    return args


# TODO: modify matchmaker

def are_there_better_motifs(n_modified_motifs, n_seqs_list, discr_exp_profile,
                            bestmi, n_bestmotif, lastmyfreq, args):

    for j in range(len(n_modified_motifs)):
        # limit it to the expressed sequences only
        # change the matchmaking procedure to incorporate degenerative nucleotides
        current_profile, time_spent = matchmaker.calculate_profile_one_motif(n_modified_motifs[j], n_seqs_list)
        myfreq = current_profile.sum() / float(len(n_seqs_list))
        tempmi = MI.mut_info(current_profile, discr_exp_profile)

        if tempmi > bestmi and current_profile.sum() > 10 and (myfreq < args.maxfreq or myfreq < lastmyfreq):
            n_bestmotif = n_modified_motifs[j].copy()
            w_bestmotif = type_conversions.n_to_w_sequence(n_bestmotif)
            bestmi = tempmi
            lastmyfreq = myfreq
            print("new motif (mi = %.4f): %s\n", bestmi, w_bestmotif.print_sequence(return_string=True))
    return bestmi, lastmyfreq, n_bestmotif


# optimize sequence of all the positions individually in random order
def optimize_motif_sequence(counter, index, n_motifs_list, n_seqs_list,
                            discr_exp_profile, active_profile,
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
                                                    n_seqs_list, discr_exp_profile,
                                                    bestmi, n_bestmotif, lastmyfreq, args)
    return bestmi, lastmyfreq, n_bestmotif


def elongate_motif(counter, index, n_motifs_list, n_seqs_list,
                    discr_exp_profile, active_profile,
                    bestmi, n_bestmotif, lastmyfreq, args):
    print("Elongating the motif %d" % counter)

    premi = bestmi

    # TODO: how to I to a first iteration of this?
    # TODO: also why premi >= bestmi, shouldn't it go the other way??

    while premi >= bestmi and n_bestmotif.length > n_motifs_list[index].length:
        n_elongated_motifs = modify_seed.elongate_motif(n_motifs_list[index])
        for j in range(len(n_elongated_motifs)):
            bestmi, lastmyfreq, n_bestmotif = are_there_better_motifs(n_elongated_motifs,
                                                          n_seqs_list, discr_exp_profile,
                                                          bestmi, n_bestmotif, lastmyfreq, args)
    return bestmi, lastmyfreq, n_bestmotif


def optimize_motifs(number_signigicant_seeds, MI_values_array, discr_exp_profile,
                    profiles_array, index_array, n_motifs_list, n_seqs_list, args, do_print = False):

    seed_indices_sorted = np.argsort(MI_values_array)[::-1]
    signif_indices = seed_indices_sorted[0 : number_signigicant_seeds]

    print("There are %d signif seeds in total" % number_signigicant_seeds)
    for counter, index in enumerate(signif_indices):
        profile = profiles_array[index]
        active_profile = profile[index_array]
        current_MI = MI_values_array[index]

        # TODO: find out if Hani's code uses parameters like doonlypositive at all

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

        n_bestmotif = n_motifs_list[index].copy()

        if do_optimize:
            # initial mi value
            init_best_mymi = MI.mut_info(active_profile, discr_exp_profile)
            lastmyfreq = active_profile.sum() / float(active_profile.shape[0])

            print("Initial MI = %.5f\n", init_best_mymi)

            bestmi, lastmyfreq, n_bestmotif = optimize_motif_sequence(counter, index,
                                    n_motifs_list, n_seqs_list,
                                    discr_exp_profile, active_profile,
                                    init_best_mymi, n_bestmotif, lastmyfreq, args)
            bestmi, lastmyfreq, n_bestmotif = elongate_motif(counter, index,
                           n_motifs_list, n_seqs_list,
                           discr_exp_profile, active_profile,
                           bestmi, n_bestmotif, lastmyfreq, args)

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

    # read occurence profiles and expression profile
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
                    n_motifs_list, n_seqs_list,
                       args, do_print=True)





if __name__ == "__main__":
    main()
