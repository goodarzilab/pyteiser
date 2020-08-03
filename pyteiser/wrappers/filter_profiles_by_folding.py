# Hani's get_dG_ratio_to_optimum filter:
# (see files folding_energy.c and matchmaker.c)
# - returns the ratio of the dG for a fold that takes the shape of the "motif" at the specified position relative to the global minimum dG for folding
# - mask motif match with Ns (function lcl_get_surrounding_sequence)
# - make an alignment matrix (Nussinov style) of the whole transcript (lcl_get_minimum_energy); calculate minimum folding energy
# - calculate folding energy for the motif (lcl_get_free_folding_energy); calculate MFE for the transcript that has the motif replaced by Ns (lcl_get_minimum_energy)
# - return the ratio: (motif energy + MFE of all the rest) / MFE of the whole transcript
# - if such ratio is smaller than a threshold value (dG_t), remove this match

# Here, I will use a similar filter. However, there will be a few distinctions
# - first of all, I don't think that MFE of the whole transcript is informative. Folding algorithms work well only in local context so I will use a window of 100 nt
# - second of all, I will be using ViennaRNA instead of dynamic programming algorithm

import os
import sys
import numpy as np
import argparse
import numba
import math
import re
import time
from subprocess import PIPE, run, Popen


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task_mapping_file", help="", type=str)

    parser.add_argument("--rna_bin_file", help="", type=str)
    parser.add_argument("--seeds_file", help="", type=str)
    parser.add_argument("--profiles_full_file", help="", type=str)
    parser.add_argument("--profiles_output_file", help="", type=str)
    parser.add_argument("--are_seeds_degenerate", help="", type=str)
    parser.add_argument("--n_profiles_check",
                        help="number of profiles to test to see if the profiles file corresponds to the seeds file",
                        type=str)
    parser.add_argument("--window_size", help="the window surrounding a match that we are folding", type=int)
    parser.add_argument("--MFE_ratio_thresh",
                        help="minimal ratio of Minimal Folding Energies for a structure with or without fixed seed match",
                        type=float)
    parser.add_argument("--indices_mode", help="compression in the index mode", type=bool)


    parser.set_defaults(
        rna_bin_file='/Users/student/Documents/hani/programs/pyteiser/data/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        seeds_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds_unique.bin',
        profiles_full_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique.bin',
        profiles_output_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique_fold_filtered.bin',

        are_seeds_degenerate='n',
        window_size=100,
        n_profiles_check=0,
        MFE_ratio_thresh=0.5,
        indices_mode=False,
    )

    args = parser.parse_args()

    return args


def compile_global_patterns():
    global energy_pattern_RNAeval
    energy_pattern_RNAeval = re.compile("\(\s*\-*\d+\.\d+\)")


def make_constraint_string(n_motif, w_motif, motif_linear_length,
                           subsequence_n, is_degenerate):
    # make one constraint string for each match
    constraint_strings = []
    # see -C parameter here: https://www.tbi.univie.ac.at/RNA/RNAfold.1.html
    subs_instances = matchmaker.find_all_motif_instances(n_motif, subsequence_n, is_degenerate=is_degenerate)
    for subs_match in subs_instances:
        structure_array = np.full(shape=subsequence_n.length,
                                  fill_value=glob_var._loop,
                                  dtype=np.uint8)
        structure_array[subs_match : subs_match + motif_linear_length] = w_motif.full_structure_encoding
        constraint_string = ''.join([glob_var._extended_structure_to_char[x] for x in structure_array])
        constraint_strings.append(constraint_string)
    return constraint_strings


def import_modules():
    package_home_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if package_home_path not in sys.path:
        sys.path.append(package_home_path)

    global structures
    global IO
    global type_conversions
    global glob_var
    global matchmaker

    import structures
    import IO
    import type_conversions
    import glob_var
    import matchmaker


def read_input_files(seeds_filename_full, rna_bin_filename):
    seqs_dict, seqs_order = IO.read_rna_bin_file(rna_bin_filename)
    w_motifs_list = IO.read_motif_file(seeds_filename_full)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]
    n_motifs_list = type_conversions.w_to_n_motifs_list(w_motifs_list)
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)

    return w_motifs_list, w_seqs_list, n_motifs_list, n_seqs_list


def make_sure_the_profile_is_correct(n_motifs_list, n_seqs_list,
                                    known_profiles_array,
                                    are_seeds_degenerate,
                                    n_to_check = 5,
                                    do_print=False):
    if are_seeds_degenerate == 'yes' or are_seeds_degenerate == 'y':
        is_degenerate = True
    else:
        is_degenerate = False

    for i in range(min(len(n_motifs_list), n_to_check)):
        motif = n_motifs_list[i]
        current_profile, time_spent = matchmaker.calculate_profile_one_motif(motif, n_seqs_list,
                                                                             is_degenerate=is_degenerate)
        if do_print:
            assert (current_profile.values == known_profiles_array[i,:]).all(), "the provided file with profiles is not matching the seeds"


def process_one_transcript_one_seed():
    pass


@numba.jit(cache=False, nopython=True, nogil=True)
def extract_subsequence(inp_seq, match_coord, motif_linear_length, window_size):
    # if a sequence is shorter than a window size, return all of it
    if inp_seq.length < window_size:
        start_coord = 0
        end_coord = inp_seq.length
    else:
        # center the window on the middle of motif; however, if motif is longer than a window, take all of it
        motif_middle = motif_linear_length // 2
        center_coord = match_coord + motif_middle
        start_coord = min((center_coord - window_size // 2), match_coord)
        end_coord = max((center_coord + window_size // 2), (center_coord + motif_middle))

        # if we are out of the transcript boundaries, shift
        if start_coord < 0:
            end_coord += (0 - start_coord)
            start_coord = 0
        elif end_coord > inp_seq.length:
            start_coord -= (end_coord - inp_seq.length)
            end_coord = inp_seq.length

    # take the subsequence
    subsequence = structures.n_sequence((end_coord - start_coord),
                                        inp_seq.nts[start_coord: end_coord])
    return subsequence


def parse_RNAfold_output(out_string):
    energy_value = 0.
    if out_string.startswith("ERROR"):
        folded_properly = False
        return (folded_properly, energy_value)
    folded_properly = True
    res = energy_pattern_RNAeval.search(out_string)
    energy_value = float(res.group()[1:-1])

    return (folded_properly, energy_value)


def compare_folding_energies(default_result, constraint_result, ratio_threshold, do_print):
    def_folded_properly, def_energy = parse_RNAfold_output(default_result)
    constr_folded_properly, constr_energy = parse_RNAfold_output(constraint_result)

    ratio_is_proper = False
    if not constr_folded_properly:
        if do_print:
            print("The sequence can't be folded with the constrains corresponding to the matching seed")
        return ratio_is_proper
    if def_energy > 0:
        print("The selected sequence does not fold properly!")
        print(default_result)
        print(def_energy)
        sys.exit(1)
    if constr_energy > 0:
        if do_print:
            print("The sequence folding with the constrains corresponding to the matching seed has positive energy")
        return ratio_is_proper
    ratio = constr_energy / def_energy
    if ratio >= ratio_threshold:
        ratio_is_proper = True
    if do_print:
        string_to_pring = "The energy ratio is %.2f. The ratio is above threshold: %r" % (ratio, ratio_is_proper)
        string_to_pring += " (energies are: %.1f, %.1f for constrained and unconstrained)" % (constr_energy, def_energy)
        print(string_to_pring)
    return ratio_is_proper


def process_one_profile_one_seed(w_motif, n_motif,
                                 w_seqs_list, n_seqs_list,
                                 profile, window_size,
                                 is_degenerate,
                                 MFE_ratio_thresh,
                                 do_print = False,
                                 do_print_subs_matches = False,
                                 do_print_progress = True,
                                 how_often_print=100):
    # prepare
    new_profile = np.zeros_like(profile)
    motif_linear_length = w_motif.linear_length
    w_motif.encode_linear_structure()
    tic = time.time()
    if do_print:
        motif_linear_sequence = w_motif.print_linear_sequence(return_string = True)
        motif_linear_structure = w_motif.print_linear_structure(return_string = True)

    # get indices of all the matched transcripts
    assert len(w_seqs_list) == profile.shape[0], "transcriptome length doesn't correspond to the profile length"
    true_indices = np.where(profile)  # get indices
    true_indices = true_indices[0]  # for some reason np.where returns a tuple
    if true_indices.shape[0] == 0:
        print("There aren't any matches in this profile")
        return profile

    # iterate through all sequences that have a match
    for k, idx in enumerate(true_indices):
        is_there_match = False
        curr_motif_instances = matchmaker.find_all_motif_instances(n_motif, n_seqs_list[idx], is_degenerate=is_degenerate)
        if do_print:
            print("Sequence %d, length %d" % (idx,n_seqs_list[idx].length))
            print("Match indices: ", ", ".join([str(x) for x in curr_motif_instances]))
        for match_coord in curr_motif_instances:
            subsequence_n = extract_subsequence(n_seqs_list[idx], match_coord, motif_linear_length, window_size)
            subsequence_w = type_conversions.n_to_w_sequence(subsequence_n)
            subsequence_string = subsequence_w.print(return_string=True)
            constraint_strings = make_constraint_string(n_motif, w_motif, motif_linear_length,
                                                       subsequence_n, is_degenerate)

            for constraint_string in constraint_strings:
                # print matches for debugging
                if do_print:
                    print("Sequence:    ", subsequence_string)
                    print("Constraints: ", constraint_string)
                if do_print and do_print_subs_matches:
                    subs_instances = matchmaker.find_all_motif_instances(n_motif, subsequence_n, is_degenerate=is_degenerate)
                    for subs_match in subs_instances:
                        print("Match (coord %d): %s" % (subs_match, subsequence_string[subs_match : subs_match + motif_linear_length]))
                        print("Motif:            %s" % motif_linear_sequence)
                        print("Structure:        %s" % motif_linear_structure)

                # fold with ViennaRNA
                # to pipe two processes, I first need to use Popen like here: https://stackoverflow.com/questions/50682514/when-python3-chain-two-subprocess-run-such-as-bash-pipe-get-error-attributeer
                # overall following the syntax from here: https://stackoverflow.com/questions/13332268/how-to-use-subprocess-command-with-pipes
                printf_default = Popen(args=["printf", "%s" % subsequence_string], stdout=PIPE, stderr=PIPE)
                RNAfold_default = run(["RNAfold"], stdin=printf_default.stdout, stdout=PIPE, stderr=PIPE, text=True)
                printf_default.terminate()
                default_result = RNAfold_default.stdout

                constrained_arguments_list = ["RNAfold", "--enforceConstraint", "-C"]
                printf_constraint = Popen(args=["printf", "%s\n%s" % (subsequence_string, constraint_string)], stdout=PIPE, stderr=PIPE)
                RNAfold_constraint = run(constrained_arguments_list, stdin=printf_constraint.stdout, stdout=PIPE, stderr=PIPE, text=True)
                printf_constraint.terminate()
                constraint_result = RNAfold_constraint.stdout

                # calculate ratio of energies and compare it to the threshold value
                ratio_is_proper = compare_folding_energies(default_result, constraint_result, MFE_ratio_thresh, do_print=do_print)
                is_there_match = is_there_match or ratio_is_proper

        if k % how_often_print == 0:
            string_to_print = "%d out of %d matches processed" % (k, true_indices.shape[0])
            if do_print_progress:
                print(string_to_print)

        new_profile[idx] = is_there_match

    toc = time.time()
    if do_print_progress:
        print("time spent on this seed: %d" % (toc-tic))

    return new_profile


def filter_profiles_by_folding(w_motifs_list, w_seqs_list,
                   n_motifs_list, n_seqs_list,
                   profiles_array,
                   output_filename,
                   window_size,
                   MFE_ratio_thresh,
                   is_degenerate,
                   do_print=False,
                   do_print_subs_matches=False,
                   do_print_progress=True,
                   how_often_print=100
                   ):
    N_seq = len(n_seqs_list)
    with open(output_filename, 'wb') as wf:
        # iterate over seeds
        for i, w_motif in enumerate(w_motifs_list):


            if i <= 5:
                continue
            if i >= 16:
                break


            n_motif = n_motifs_list[i]
            current_profile = profiles_array[i, :]
            filtered_profile = process_one_profile_one_seed(w_motif, n_motif,
                                         w_seqs_list, n_seqs_list,
                                         current_profile, window_size, is_degenerate,
                                         MFE_ratio_thresh,
                                         do_print=do_print,
                                         do_print_subs_matches=do_print_subs_matches,
                                         do_print_progress=do_print_progress,
                                         how_often_print=how_often_print)

            filtered_profile_w = structures.w_profile(N_seq)
            filtered_profile_w.values = filtered_profile
            filtered_profile_w.compress_indices()
            wf.write(filtered_profile_w.bytestring_indices)

            if do_print:
                difference_u_f = np.logical_and(current_profile, np.invert(filtered_profile))
                print("%d out of %d transcripts were filtered out: " % (difference_u_f.sum(), current_profile.sum()))



def main():
    import_modules()
    compile_global_patterns()
    args = handler()

    w_motifs_list, w_seqs_list, n_motifs_list, n_seqs_list  = read_input_files(args.seeds_file, args.rna_bin_file)
    decompressed_profiles_array = IO.unpack_profiles_file(args.profiles_full_file,
                                                          args.indices_mode,
                                                          do_print=True)

    if args.n_profiles_check > 0:
        make_sure_the_profile_is_correct(n_motifs_list, n_seqs_list,
                                         decompressed_profiles_array,
                                         args.are_seeds_degenerate,
                                         n_to_check = args.n_profiles_check,
                                         do_print=True)

    if args.are_seeds_degenerate == 'yes' or args.are_seeds_degenerate == 'y':
        are_seeds_degenerate = True
    else:
        are_seeds_degenerate = False

    filter_profiles_by_folding(w_motifs_list, w_seqs_list,
                   n_motifs_list, n_seqs_list,
                   decompressed_profiles_array,
                   args.profiles_output_file,
                   window_size = args.window_size,
                   MFE_ratio_thresh = args.MFE_ratio_thresh,
                   is_degenerate = are_seeds_degenerate)



if __name__ == "__main__":
    main()
