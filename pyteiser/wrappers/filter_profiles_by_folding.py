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
# - second of all, I will try using ViennaRNA; if it turns out to be too slow, I will use a dynamic algorithm as well

import os
import sys
import numpy as np
import argparse
import numba
import math
import subprocess


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task_mapping_file", help="", type=str)

    parser.add_argument("--rna_bin_file", help="", type=str)
    parser.add_argument("--seeds_file", help="", type=str)
    parser.add_argument("--profiles_full_file", help="", type=str)
    parser.add_argument("--are_seeds_degenerate", help="", type=str)
    parser.add_argument("--n_profiles_check",
                        help="number of profiles to test to see if the profiles file corresponds to the seeds file",
                        type=str)
    parser.add_argument("--window_size", help="the window surrounding a match that we are folding", type=int)

    parser.set_defaults(
        rna_bin_file='/Users/student/Documents/hani/programs/pyteiser/data/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        seeds_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds_unique.bin',
        profiles_full_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique.bin',

        are_seeds_degenerate='n',
        window_size=100,
        n_profiles_check=0,

    )

    args = parser.parse_args()

    return args


def make_constraint_string(n_motif, w_motif, motif_linear_length,
                           subsequence_n, is_degenerate):
    # see -C parameter here: https://www.tbi.univie.ac.at/RNA/RNAfold.1.html
    structure_array = np.full(shape=subsequence_n.length,
                               fill_value=glob_var._loop,
                               dtype=np.uint8)

    subs_instances = matchmaker.find_all_motif_instances(n_motif, subsequence_n, is_degenerate=is_degenerate)
    for subs_match in subs_instances:
        structure_array[subs_match : subs_match + motif_linear_length] = w_motif.full_structure_encoding
    constraint_string = ''.join([glob_var._extended_structure_to_char[x] for x in structure_array])
    return constraint_string


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


def process_one_profile_one_seed(w_motif, n_motif,
                                 w_seqs_list, n_seqs_list,
                                 profile, window_size,
                                 is_degenerate,
                                 do_print = True,
                                 do_print_subs_matches = True):
    # prepare
    motif_linear_length = w_motif.linear_length
    w_motif.encode_linear_structure()
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
        curr_motif_instances = matchmaker.find_all_motif_instances(n_motif, n_seqs_list[idx], is_degenerate=is_degenerate)
        if do_print:
            print("Sequence %d, length %d" % (idx,n_seqs_list[idx].length))
            print("Match indices: ", ", ".join([str(x) for x in curr_motif_instances]))
        for match_coord in curr_motif_instances:
            subsequence_n = extract_subsequence(n_seqs_list[idx], match_coord, motif_linear_length, window_size)
            subsequence_w = type_conversions.n_to_w_sequence(subsequence_n)
            subsequence_string = subsequence_w.print(return_string=True)
            constraint_string = make_constraint_string(n_motif, w_motif, motif_linear_length,
                                                       subsequence_n, is_degenerate)
            if do_print:
                print("Sequence:    ", subsequence_string)
                print("Constraints: ", constraint_string)
            if do_print and do_print_subs_matches:
                subs_instances = matchmaker.find_all_motif_instances(n_motif, subsequence_n, is_degenerate=is_degenerate)
                for subs_match in subs_instances:
                    print("Match (coord %d): %s" % (subs_match, subsequence_string[subs_match : subs_match + motif_linear_length]))
                    print("Motif:            %s" % motif_linear_sequence)
                    print("Structure:        %s" % motif_linear_structure)

            viennarna_command = ""




        # analog of get_dG_ratio_to_optimum

        if k == 10:
            break


    # this function only works with n_motif and n_sequence classes,





def fold_instances(w_motifs_list, w_seqs_list,
                   n_motifs_list, n_seqs_list,
                   profiles_array, window_size,
                   is_degenerate):
    # iterate over instances
    #true_indices = np.where(self.values)  # get indices
    #true_indices = true_indices[0]  # for some reason np.where returns a tuple

    for i, w_motif in enumerate(w_motifs_list):
        n_motif = n_motifs_list[i]
        current_profile = profiles_array[i, :]
        process_one_profile_one_seed(w_motif, n_motif, w_seqs_list, n_seqs_list,
                                     current_profile, window_size, is_degenerate)
        break


def main():
    import_modules()
    args = handler()

    w_motifs_list, w_seqs_list, n_motifs_list, n_seqs_list  = read_input_files(args.seeds_file, args.rna_bin_file)
    decompressed_profiles_array = IO.unpack_profiles_file(args.profiles_full_file, do_print=True)

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

    fold_instances(w_motifs_list, w_seqs_list,
                   n_motifs_list, n_seqs_list,
                   decompressed_profiles_array,
                   window_size = args.window_size,
                   is_degenerate = are_seeds_degenerate)



if __name__ == "__main__":
    main()
