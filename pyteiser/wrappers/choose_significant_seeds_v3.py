import numpy as np
import argparse
import math

import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

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
    parser.add_argument("--task_mapping_file", help="", type=str)

    parser.add_argument("--profiles_folder", help="", type=str)
    parser.add_argument("--MI_values_folder", help="", type=str)
    parser.add_argument("--passed_seed_folder", help="", type=str)
    parser.add_argument("--passed_profiles_folder", help="", type=str)
    parser.add_argument("--seed_folder", help="", type=str)

    parser.add_argument("--profiles_filename_template", help="", type=str)
    parser.add_argument("--MI_values_filename_template", help="", type=str)
    parser.add_argument("--passed_seed_filename_template", help="", type=str)
    parser.add_argument("--passed_profiles_filename_template", help="", type=str)
    parser.add_argument("--seed_filename_template", help="", type=str)

    parser.add_argument("--rna_bin_file", help="", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)

    parser.add_argument("--n_permutations", help="number of permutations for the rnak test for a seed", type=int)
    parser.add_argument("--max_pvalue", help="maximal acceptable p-value", type=float)
    parser.add_argument("--min_zscore", help="maximal acceptable p-value", type=float)
    parser.add_argument("--step_1_jump", help="", type=int)
    parser.add_argument("--step_2_min_interval", help="", type=int)
    parser.add_argument("--step_1_min_fraction", help="", type=float)
    parser.add_argument("--step_2_min_fraction", help="", type=float)
    parser.add_argument("--step_3_min_fraction", help="", type=float)

    parser.add_argument("--print_qstat", help="", type=str)
    parser.add_argument("--path_to_qstat", help="", type=str)
    parser.add_argument("--indices_mode", help="compression in the index mode", type=bool)
    parser.add_argument("--index_bit_width", help="number of bits per one index when compressing", type=bool)



    parser.set_defaults(
        seed_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/seeds/seeds_4-7_4-9_4-6_14-20/motifs_per_file_30k',
        profiles_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/profiles/profiles_4-7_4-9_4-6_14-20/profiles_per_file_30k',
        MI_values_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/MI_values/MI_values_4-7_4-9_4-6_14-20/MI_values_per_file_30k',
        passed_seed_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_seed/passed_seed_4-7_4-9_4-6_14-20_individual/passed_seed_30k',
        passed_profiles_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_individual/passed_profiles_30k',
        seed_filename_template='seeds_4-7_4-9_4-6_14-20_30k',
        profiles_filename_template='profiles_4-7_4-9_4-6_14-20_30k',
        MI_values_filename_template='MI_values_4-7_4-9_4-6_14-20_30k',
        passed_seed_filename_template='passed_seed_4-7_4-9_4-6_14-20_30k',
        passed_profiles_filename_template='passed_profiles_4-7_4-9_4-6_14-20_30k',

        rna_bin_file='/wynton/home/goodarzi/khorms/pyteiser_root/data/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        exp_mask_file='/wynton/home/goodarzi/khorms/pyteiser_root/data/mask_files/TARBP2_decay_t_score_mask.bin',

        n_permutations = 1000, # takes 1 second per 100 permutations, Hani's default number of permutations is 1*10^6
        max_pvalue = 0.0001, # Hani's default threshold is 1*10^-7
        min_zscore = -1,

        step_1_jump = 100, # Hani's default jump is 200
        step_2_min_interval=10,

        step_1_min_fraction = 0.8,
        step_2_min_fraction= 0.8,
        step_3_min_fraction= 0.9,

        path_to_qstat='/opt/sge/bin/lx-amd64/qstat',
        print_qstat='y',
        indices_mode=False,
        index_bit_width=24,

    )

    args = parser.parse_args()

    return args


def import_modules():
    current_wd = os.getenv('SGE_O_WORKDIR')
    subpackage_folder_path = os.path.abspath(os.path.join(current_wd, '..'))
    if subpackage_folder_path not in sys.path:
        sys.path.append(subpackage_folder_path)

    global structures
    global MI
    global IO
    global sge
    global statistic_tests

    import structures
    import MI
    import IO
    import sge
    import statistic_tests


def get_current_in_out_filenames(args, env_variables_dict, mapping_dict):
    file_index_to_use =  mapping_dict[env_variables_dict["task_id"]]

    profiles_filename_short = "%s_%s.bin" % (args.profiles_filename_template, file_index_to_use)
    MI_values_filename_short = "%s_%s.bin" % (args.MI_values_filename_template, file_index_to_use)
    passed_seed_filename_short = "%s_%s.bin" % (args.passed_seed_filename_template, file_index_to_use)
    passed_profiles_filename_short = "%s_%s.bin" % (args.passed_profiles_filename_template, file_index_to_use)
    seed_filename_short = "%s_%s.bin" % (args.seed_filename_template, file_index_to_use)

    profiles_filename_full = os.path.join(args.profiles_folder, profiles_filename_short)
    MI_values_filename_full = os.path.join(args.MI_values_folder, MI_values_filename_short)
    passed_seed_filename_full = os.path.join(args.passed_seed_folder, passed_seed_filename_short)
    passed_profiles_filename_full = os.path.join(args.passed_profiles_folder, passed_profiles_filename_short)
    seed_filename_full = os.path.join(args.seed_folder, seed_filename_short)

    rna_bin_filename = args.rna_bin_file
    exp_mask_filename = args.exp_mask_file

    return profiles_filename_full, MI_values_filename_full, \
           passed_seed_filename_full, passed_profiles_filename_full, \
           seed_filename_full, \
           rna_bin_filename, exp_mask_filename


def find_fraction_denominator(x):
    rounded_denominator = round(1 / (1 - x), 4)
    return math.ceil(rounded_denominator)


def get_current_statistics(index, MI_values_array, profiles_array,
                           index_array, discr_exp_profile, nbins, args):
    profile = profiles_array[index]
    active_profile = profile[index_array]
    current_MI = MI_values_array[index]

    if current_MI == -1:
        return args.max_pvalue + 0.1, args.min_zscore - 0.1

    assert (np.isclose(current_MI, MI.mut_info(active_profile, discr_exp_profile, x_bins=2, y_bins=nbins), rtol=1e-10))

    pvalue, z_score = statistic_tests.MI_get_pvalue_and_zscore(active_profile, discr_exp_profile, nbins,
                                                               current_MI, args.n_permutations)
    return pvalue, z_score


def check_one_seed(index, counter, seed_pass, MI_values_array,
                   profiles_array, index_array, discr_exp_profile, nbins,
                   args, do_print):
    pv_defined = False
    if seed_pass[index] != 0:  # if we have checked this seed before
        check = seed_pass[index]
    else:
        pvalue, z_score = get_current_statistics(index, MI_values_array, profiles_array,
                                                 index_array, discr_exp_profile, nbins, args)
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
        discr_exp_profile, nbins, profiles_array, index_array,
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
                                       profiles_array, index_array, discr_exp_profile, nbins,
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
                                 discr_exp_profile, nbins, profiles_array, index_array,
                                args, do_print = False):
    first_negative_seed = -1

    counter = 0
    is_over_threshold = False

    while (counter < len(seed_indices_sorted)) and not is_over_threshold:
        is_over_threshold, seed_pass = check_N_consecutive(args.step_1_min_fraction,
                            counter, MI_values_array, seed_indices_sorted, seed_pass,
                            discr_exp_profile, nbins, profiles_array, index_array,
                            args, do_print)
        first_negative_seed = counter
        counter += args.step_1_jump

    if not is_over_threshold:
        print("No lower boundary for threshold found!")
        return len(MI_values_array) - 1, seed_pass

    return first_negative_seed, seed_pass


def step_2_decreasing_intervals(first_negative_seed, MI_values_array, seed_indices_sorted,
                         profiles_array, index_array, discr_exp_profile, nbins, seed_pass,
                         do_print, args):
    last_positive_seed = first_negative_seed - args.step_1_jump # this is where the last positive seed was

    upper_boundary = first_negative_seed - args.step_1_jump # upper limit: first negative seed - one jump (just in case)
    upper_boundary = min(upper_boundary, 0) # upper boundary shouldn't be less than zero
    lower_boundary = min(len(MI_values_array) - 1,
                        first_negative_seed + args.step_1_jump) # lower limit
                            # the first negative seed + step_1_jump
                            # or, if there were none, the end of the profiles list



    while (lower_boundary - upper_boundary) > args.step_2_min_interval:
        counter = upper_boundary + (lower_boundary - upper_boundary) // 2
        is_over_threshold, seed_pass = check_N_consecutive(args.step_1_min_fraction,
                                                           counter, MI_values_array, seed_indices_sorted, seed_pass,
                                                           discr_exp_profile, nbins, profiles_array, index_array,
                                                           args, do_print)
        if not is_over_threshold:
            # some seed passed, go down half interval
            upper_boundary = counter
            last_positive_seed = counter
        else:
            # enough seeds didn't pass, go up half interval
            lower_boundary = counter

    last_positive_seed = max(last_positive_seed, -1) # last positive seed could become negative if the seed 0 didn't pass
                                                     # at the first stage
    return last_positive_seed, seed_pass


def get_current_fraction_from_array(x):
    number_not_passed = np.logical_not(x).sum()
    return float(number_not_passed) / x.shape[0]


def step_3_confirm_consec_not_passing_seeds(last_positive_seed, MI_values_array, seed_indices_sorted,
                         profiles_array, index_array, discr_exp_profile, nbins, seed_pass,
                         do_print, args):
    denominator = find_fraction_denominator(args.step_3_min_fraction)
    if last_positive_seed + denominator > len(MI_values_array):
        print("Reached the end of file, no threshold found! Suppose that all the seeds are significant")
        sys.exit(1)

    current_passed_array = np.zeros(denominator, dtype=np.bool)
    for i in range(denominator):
        counter = i + last_positive_seed + 1
        index = seed_indices_sorted[counter]
        did_seed_pass, seed_pass = check_one_seed(index, counter, seed_pass, MI_values_array,
                                       profiles_array, index_array, discr_exp_profile, nbins,
                                       args, do_print)
        current_passed_array[i] = did_seed_pass

    current_fraction = get_current_fraction_from_array(current_passed_array)
    fraction_exceeded = (current_fraction > args.step_3_min_fraction)

    counter = last_positive_seed + denominator

    while counter < len(MI_values_array) and not fraction_exceeded:
        index = seed_indices_sorted[counter]
        did_seed_pass, seed_pass = check_one_seed(index, counter, seed_pass, MI_values_array,
                                                  profiles_array, index_array, discr_exp_profile, nbins,
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


def determine_mi_threshold(MI_values_array, discr_exp_profile, nbins,
                           profiles_array, index_array, args, do_print = False):

    seed_indices_sorted = np.argsort(MI_values_array)[::-1]

    seed_pass = np.zeros(MI_values_array.shape[0], dtype=np.int8) # zero means no info

    if do_print:
        print("Find the lower boundary for the threshold")
    first_negative_seed, seed_pass = step_1_determine_thresh_lower_limit(MI_values_array, seed_indices_sorted, seed_pass,
                                                     discr_exp_profile, nbins, profiles_array, index_array,
                                                    args, do_print)
    if do_print:
        print("The first seed that didn't pass is: ", first_negative_seed, '\n')
        print("Decreasing intervals phase")
    last_positive_seed, seed_pass = step_2_decreasing_intervals(first_negative_seed, MI_values_array, seed_indices_sorted,
                                                     profiles_array, index_array, discr_exp_profile, nbins,
                                                     seed_pass, do_print, args)
    if do_print:
        print("The last seed that passed is: ", last_positive_seed, '\n')
        print("Find 10 consecutive seeds that don't pass")

    last_positive_seed, seed_pass = step_3_confirm_consec_not_passing_seeds(last_positive_seed, MI_values_array,
                                            seed_indices_sorted, profiles_array, index_array, discr_exp_profile, nbins,
                                            seed_pass, do_print, args)

    if do_print:
        print("The last seed that passed is: ", last_positive_seed, '\n')

    return last_positive_seed

    # IO.write_seed_significancy_threshold(last_positive_seed, threshold_filename)


def write_seeds_passed(last_positive_seed, MI_values_array, w_motifs_list,
                       passed_seed_filename):
    if last_positive_seed < 0:
        total_bitstring = np.uint32(0).tobytes()
    else:
        seed_indices_sorted = np.argsort(MI_values_array)[::-1]
        indices_passed = seed_indices_sorted[0 : last_positive_seed + 1]
        seeds_passed_list = [w_motifs_list[x] for x in indices_passed]

        seeds_bitstrings = []

        for motif in seeds_passed_list:
            motif.compress()
            seeds_bitstrings.append(motif.bytestring)

        total_bitstring = np.uint32(len(seeds_bitstrings)).tobytes() + b''.join(seeds_bitstrings)

    with open(passed_seed_filename, 'wb') as wf:
        wf.write(total_bitstring)


def write_profiles_passed(last_positive_seed, MI_values_array, profiles_array,
                        passed_profiles_filename, indices_mode, index_bit_width):
    if last_positive_seed < 0:
        total_bitstring = np.uint32(0).tobytes()
    else:
        seed_indices_sorted = np.argsort(MI_values_array)[::-1]
        indices_passed = seed_indices_sorted[0 : last_positive_seed + 1]

        profiles_passed_list = profiles_array[indices_passed]
        profiles_bitstrings = []

        for i in range(profiles_passed_list.shape[0]):
            current_profile = structures.w_profile(profiles_passed_list[i].shape[0])
            current_profile.values = profiles_passed_list[i]

            if indices_mode:
                current_profile.compress_indices(width = index_bit_width)
                profiles_bitstrings.append(current_profile.bytestring_indices)
            else:
                current_profile.compress()
                profiles_bitstrings.append(current_profile.bytestring)

        total_bitstring = np.uint32(len(profiles_bitstrings)).tobytes() + b''.join(profiles_bitstrings)

    with open(passed_profiles_filename, 'wb') as wf:
        wf.write(total_bitstring)


def main():
    # I only import things if I run this script itself
    # do relative import based on current working directory
    # otherwise I have to install the package for relative import to work
    import_modules()

    args = handler()

    # get mapping of task ids to input files
    mapping_dict = sge.parse_task_mapping_file(args.task_mapping_file)
    # get the task id
    env_variables_dict = sge.get_env_variables()

    # get the names of input and output files
    profiles_filename_full, MI_values_filename_full, \
    passed_seed_filename_full, passed_profiles_filename, \
    seed_filename_full, \
    rna_bin_filename, exp_mask_filename = get_current_in_out_filenames(args, env_variables_dict, mapping_dict)

    # read motifs, their profiles and MI values
    profiles_array, index_array, values_array = IO.unpack_profiles_and_mask(profiles_filename_full,
                                                                         exp_mask_filename,
                                                                         args.indices_mode,
                                                                         do_print=True)
    w_motifs_list = IO.read_motif_file(seed_filename_full)
    MI_values_array, nbins = IO.read_MI_values(MI_values_filename_full)

    # find the threshold
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins)
    last_positive_seed = determine_mi_threshold(MI_values_array, discr_exp_profile, nbins,
                           profiles_array, index_array, args, do_print = True)

    write_seeds_passed(last_positive_seed, MI_values_array, w_motifs_list,
                       passed_seed_filename_full)
    write_profiles_passed(last_positive_seed, MI_values_array, profiles_array,
                          passed_profiles_filename, args.indices_mode, args.index_bit_width)


    if args.print_qstat == 'y':
        sge.print_qstat_proc(env_variables_dict, args.path_to_qstat)



if __name__ == "__main__":
    main()
