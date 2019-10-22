import numpy as np
import argparse
import os
import sys
import math
import copy


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--unique_seeds_filename", help="best representatives of each family", type=str)
    parser.add_argument("--unique_profiles_filename", help="profiles of best representatives of each family",
                                                                                             type=str)

    parser.add_argument("--rna_bin_file", help="referense transcriptome in binary format", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)


    parser.add_argument("--optimized_seeds_folder", help="output: optimized seeds", type=str)
    parser.add_argument("--optimized_profiles_folder", help="output: profiles of optimized seeds", type=str)
    parser.add_argument("--optimized_MI_pv_zscores_folder", help="output: MI values, p-values and z-scores", type=str)
    parser.add_argument("--robustness_array_folder", help="output: vector indicating which seeds have passed the robustness test", type=str)

    parser.add_argument("--optimized_seeds_filename_template", help="", type=str)
    parser.add_argument("--optimized_profiles_filename_template", help="", type=str)
    parser.add_argument("--optimized_MI_pv_zscores_filename_template", help="", type=str)
    parser.add_argument("--robustness_array_filename_template", help="", type=str)


    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)
    parser.add_argument("--maxfreq", help="maximal seed frequency in the sequences analyzed", type=float)
    parser.add_argument("--min_occurences", help="minimal number of seed occurence in the transcriptome"
                                                 " for a seed to be considered", type=int)
    parser.add_argument("--random_noseed", help="when choosing the order of positions to optimize, "
                                                "do not set the random number generator to a specific seed", type=bool)
    parser.add_argument("--jackknife_n_permutations", help="number of permutations for pvalue calculation in "
                                                           "jackknife test", type=int)
    parser.add_argument("--jackknife_max_pvalue", help="maximal pvalue for jackknife test", type=float)
    parser.add_argument("--jackknife_n_samples", help="how many permutations to do in jackknife test", type=int)
    parser.add_argument("--jackknife_fraction_retain", help="what fraction of the sample to retain for each test",
                                                                                    type=float)
    parser.add_argument("--jackknife_min_fraction_passed", help="what fraction of all iterations should"
                                                                "pass to consider the motif robust", type=float)

    parser.add_argument("--size_of_chunks", help="how many seeds should 1 process take on", type=float)

    parser.set_defaults(
        # unique_seeds_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds_unique.bin",
        # unique_profiles_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique.bin",
        #
        # optimized_seeds_folder='/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined',
        # optimized_profiles_folder='/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined',
        # optimized_MI_pv_zscores_folder='/Users/student/Documents/hani/programs/pyteiser/data/optimized_seeds_characteristics/seeds_4-7_4-9_4-6_14-20_individual',
        # robustness_array_folder='/Users/student/Documents/hani/programs/pyteiser/data/seeds_robustness/seeds_4-7_4-9_4-6_14-20_individual',
        #
        # optimized_seeds_filename_template='test_1_2_seeds_optimized',
        # optimized_profiles_filename_template='test_1_2_profiles_optimized',
        # optimized_MI_pv_zscores_filename_template='test_1_2_characteristics',
        # robustness_array_filename_template='test_1_2_robustness',
        #
        # rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        # exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',

        nbins=15,
        maxfreq = 0.5, # default value from Hani's program is 0.5
        min_occurences=10,
        n_permutations=1000,  # takes 1 second per 100 permutations, Hani's default number of permutations is 1*10^6
        random_noseed=0,
        jackknife_n_samples = 10,
        jackknife_fraction_retain = 0.66,
        jackknife_n_permutations=1000,
        jackknife_max_pvalue=0.0001,
        jackknife_min_fraction_passed = 0.6,

        size_of_chunks=10,
    )

    args = parser.parse_args()

    return args


def import_modules():
    current_wd = os.getenv('SGE_O_WORKDIR')
    subpackage_folder_path = os.path.abspath(os.path.join(current_wd, '..'))

    if subpackage_folder_path not in sys.path:
        sys.path.append(subpackage_folder_path)

    global MI
    global IO
    global sge
    global structures
    global modify_seed
    global type_conversions
    global matchmaker
    global statistic_tests

    import MI
    import IO
    import sge
    import structures
    import modify_seed
    import type_conversions
    import matchmaker
    import statistic_tests



def chunk_up_input_files(seeds_initial, profiles_initial, size_of_chunks):
    seeds_number = len(seeds_initial)
    print("Starting with %d initial seeds" % seeds_number)

    number_of_chunks = math.ceil(seeds_number / size_of_chunks)

    profiles_chunks = np.array_split(profiles_initial, number_of_chunks)
    seed_array = np.array(seeds_initial)
    seed_chunks = np.array_split(seed_array, number_of_chunks)

    for pr_ch, s_ch in zip(profiles_chunks, seed_chunks):
        assert(len(pr_ch) == len(s_ch))

    return seed_chunks, profiles_chunks


def pick_one_chunk(seed_chunks, profiles_chunks, env_variables_dict):
    chunk_number = int(env_variables_dict["task_id"]) - 1
    print("Processing the chunk number ", chunk_number)
    seed_right_chunk = seed_chunks[chunk_number]
    profiles_right_chunk = profiles_chunks[chunk_number]

    return seed_right_chunk, profiles_right_chunk


def are_there_better_motifs(n_modified_motifs, seqs_of_interest, discr_exp_profile, nbins,
                            bestmi, n_bestmotif, lastmyfreq, args, do_print = True):

    for curr_motif in n_modified_motifs:
        current_profile, time_spent = matchmaker.calculate_profile_one_motif(curr_motif,
                                                                             seqs_of_interest,
                                                                            is_degenerate = True)
        myfreq = current_profile.values.sum() / float(len(seqs_of_interest))
        tempmi = MI.mut_info(current_profile.values, discr_exp_profile, x_bins=2, y_bins=nbins)

        if tempmi > bestmi and current_profile.sum() > args.min_occurences and (myfreq < args.maxfreq or myfreq < lastmyfreq):
            n_bestmotif = structures.copy_n_motif(curr_motif)
            w_bestmotif = type_conversions.n_to_w_motif(n_bestmotif)
            bestmi = tempmi
            lastmyfreq = myfreq
            if do_print:
                print("New motif (MI = %.4f): %s" % (bestmi, w_bestmotif.print_sequence(return_string=True)))
                # w_bestmotif.print()
                # w_bestmotif.print_linear()
                #print("Current frequency: %.4f" % lastmyfreq)
    return bestmi, lastmyfreq, n_bestmotif


# optimize sequence of all the positions individually in random order
def optimize_motif_sequence(n_bestmotif, init_best_MI, seqs_of_interest,
                            discr_exp_profile, nbins, lastmyfreq, args,
                            do_print = False, random_noseed=False):
    bestmi = init_best_MI

    if random_noseed:
        np.random.seed(1543)

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


def elongate_motif(n_bestmotif, init_best_MI, seqs_of_interest,
                   discr_exp_profile, nbins, lastmyfreq,
                   args, do_print = False):
    bestmi = init_best_MI

    keep_elongating = True

    # simple emulator of do {} while {} in python: https://stackoverflow.com/questions/743164/emulate-a-do-while-loop-in-python
    while keep_elongating:
        n_elongated_motifs = modify_seed.elongate_motif(n_bestmotif)

        old_best_mi = bestmi
        old_best_motif = structures.copy_n_motif(n_bestmotif)

        new_bestmi, lastmyfreq, n_bestmotif = are_there_better_motifs(n_elongated_motifs,
                                                    seqs_of_interest, discr_exp_profile, nbins,
                                                    bestmi, n_bestmotif, lastmyfreq, args,
                                                    do_print = do_print)

        keep_elongating = ((new_bestmi >= old_best_mi) and
                           (n_bestmotif.length > old_best_motif.length))

    return bestmi, lastmyfreq, n_bestmotif


def get_characteristics(n_bestmotif, seqs_of_interest,
                        discr_exp_profile, nbins, args,
                        do_print = False):
    bestmotif_profile, _time = matchmaker.calculate_profile_one_motif(n_bestmotif, seqs_of_interest,
                                                                      is_degenerate = True)
    bestmotif_mi = MI.mut_info(bestmotif_profile.values, discr_exp_profile, x_bins=2, y_bins=nbins)
    pvalue, z_score = statistic_tests.MI_get_pvalue_and_zscore(bestmotif_profile.values, discr_exp_profile, nbins,
                                                               bestmotif_mi, args.n_permutations)
    if do_print:
        print("The final p-value is: %.4f, z-score is: %.3f" % (pvalue, z_score))
    return bestmotif_profile, bestmotif_mi, pvalue, z_score


def check_robustness(bestmotif_profile,
                    discr_exp_profile, nbins, args,
                    do_print = False):


    passed_jacknife = statistic_tests.jackknife_test(
                        bestmotif_profile.values, discr_exp_profile, nbins,
                        args.jackknife_n_permutations,
                        args.jackknife_max_pvalue,
                        args.jackknife_n_samples,
                        args.jackknife_fraction_retain,
                        args.jackknife_min_fraction_passed,
                        do_print = do_print)


    return passed_jacknife



def optimize_motifs(seeds_initial, profiles_initial,
                    discr_exp_profile, nbins, index_array, seqs_of_interest,
                    args, do_print = True):
    seeds_optimized = copy.deepcopy(seeds_initial)
    profiles_optimized = np.zeros((len(seeds_initial), discr_exp_profile.shape[0]), dtype=bool)
    # seed_charact_array keeps MI values, p-values and z-scores
    seed_charact_array = np.zeros((len(seeds_initial), 3), dtype=np.float64)
    robustness_array = np.zeros(len(seeds_initial), dtype=bool)

    for i, motif in enumerate(seeds_initial):
        profile = profiles_initial[i]
        active_profile = profile[index_array]
        n_bestmotif = type_conversions.w_to_n_motif(seeds_initial[i])

        # initial mi value
        init_best_MI = MI.mut_info(active_profile, discr_exp_profile, x_bins=2, y_bins=nbins)
        lastmyfreq = active_profile.sum() / float(active_profile.shape[0])

        if do_print:
            w_bestmotif = type_conversions.n_to_w_motif(n_bestmotif)
            print("Optimzing the sequence of motif %d (sequence is %s). Initial MI = %.5f" %
                            (i, w_bestmotif.print_sequence(return_string=True), init_best_MI))
            #print("Initial frequency: %.4f" % lastmyfreq)

        bestmi, lastmyfreq, n_bestmotif = optimize_motif_sequence(n_bestmotif, init_best_MI, seqs_of_interest,
                            discr_exp_profile, nbins, lastmyfreq, args, do_print = do_print,
                            random_noseed = args.random_noseed)

        if do_print:
            print("Elongating motif %d" % i)

        bestmi, lastmyfreq, n_bestmotif = elongate_motif(n_bestmotif, bestmi, seqs_of_interest,
                            discr_exp_profile, nbins, lastmyfreq, args, do_print = do_print)

        w_bestmotif = type_conversions.n_to_w_motif(n_bestmotif)
        bestmotif_profile, bestmotif_mi, pvalue, z_score = get_characteristics(
                                                            n_bestmotif, seqs_of_interest,
                                                            discr_exp_profile, nbins, args,
                                                            do_print=do_print)

        if do_print:
            print("Checking robustness of the optimized motif %d (sequence %s)" %
                  (i, w_bestmotif.print_sequence(return_string=True)))

        is_robust = check_robustness(bestmotif_profile,
                                    discr_exp_profile, nbins, args,
                                    do_print = do_print)

        seeds_optimized[i] = w_bestmotif
        profiles_optimized[i] = bestmotif_profile.values
        seed_charact_array[i, : ] = np.array([bestmotif_mi, pvalue, z_score], dtype=np.float64)
        robustness_array[i] = is_robust

    return seeds_optimized, profiles_optimized, \
           seed_charact_array, robustness_array


def read_sequences(rna_bin_filename):
    seqs_dict, seqs_order = IO.read_rna_bin_file(rna_bin_filename)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)

    return n_seqs_list


def make_output_filenames(env_variables_dict, args):
    file_index_to_use =  env_variables_dict["task_id"]

    seed_filename_short = "%s_%s.bin" % (args.optimized_seeds_filename_template, file_index_to_use)
    profiles_filename_short = "%s_%s.bin" % (args.optimized_profiles_filename_template, file_index_to_use)
    char_filename_short = "%s_%s.bin" % (args.optimized_MI_pv_zscores_filename_template, file_index_to_use)
    robustness_filename_short = "%s_%s.bin" % (args.robustness_array_filename_template, file_index_to_use)

    seeds_filename_full = os.path.join(args.optimized_seeds_folder, seed_filename_short)
    profiles_filename_full = os.path.join(args.optimized_profiles_folder, profiles_filename_short)
    char_filename_full = os.path.join(args.optimized_MI_pv_zscores_folder, char_filename_short)
    robustness_filename_full = os.path.join(args.robustness_array_folder, robustness_filename_short)

    return seeds_filename_full, profiles_filename_full, \
           char_filename_full, robustness_filename_full



def main():
    import_modules()
    args = handler()

    n_seqs_list = read_sequences(args.rna_bin_file)
    index_array, values_array = IO.unpack_mask_file(args.exp_mask_file)
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins = args.nbins)
    seeds_initial = IO.read_motif_file(args.unique_seeds_filename)
    profiles_initial = IO.unpack_profiles_file(args.unique_profiles_filename)
    seqs_of_interest = [n_seqs_list[x] for x in range(index_array.shape[0]) if index_array[x]]

    # get the task id
    env_variables_dict = sge.get_env_variables()
    seed_chunks, profiles_chunks = chunk_up_input_files(seeds_initial, profiles_initial, args.size_of_chunks)
    seed_right_chunk, profiles_right_chunk = pick_one_chunk(seed_chunks, profiles_chunks, env_variables_dict)

    seeds_filename_full, profiles_filename_full, \
    char_filename_full, robustness_filename_full = make_output_filenames(env_variables_dict, args)


    seeds_optimized, profiles_optimized, \
    seed_charact_array, robustness_array  = optimize_motifs(seed_right_chunk, profiles_right_chunk,
                                            discr_exp_profile, args.nbins, index_array, seqs_of_interest,
                                            args, do_print=True)

    IO.write_list_of_seeds(seeds_optimized, seeds_filename_full)
    IO.write_array_of_profiles(profiles_optimized, profiles_filename_full)
    IO.write_np_array(seed_charact_array, char_filename_full)
    IO.write_np_array(robustness_array, robustness_filename_full)




if __name__ == "__main__":
    main()
