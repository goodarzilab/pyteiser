import os
import sys
import numpy as np
import argparse
import math
import subprocess


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task_mapping_file", help="", type=str)

    parser.add_argument("--rna_bin_file", help="", type=str)
    parser.add_argument("--seeds_file", help="", type=str)
    parser.add_argument("--profiles_full_file", help="", type=str)
    parser.add_argument("--are_seeds_degenerate", help="", type=str)
    parser.add_argument("--n_profiles_check", help="number of profiles to test to see if the profiles file corresponds to the seeds file", type=str)

    parser.set_defaults(
        rna_bin_file='/Users/student/Documents/hani/programs/pyteiser/data/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        seeds_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds_unique.bin',
        profiles_full_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique.bin',

        are_seeds_degenerate='n',
        n_profiles_check=1,

    )

    args = parser.parse_args()

    return args


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

    return n_motifs_list, n_seqs_list


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


def fold_instances():
    # iterate over instances
    true_indices = np.where(self.values)  # get indices
    true_indices = true_indices[0]  # for some reason np.where returns a tuple


def main():
    import_modules()
    args = handler()

    n_motifs_list, n_seqs_list = read_input_files(args.seeds_file, args.rna_bin_file)
    decompressed_profiles_array = IO.unpack_profiles_file(args.profiles_full_file, do_print=True)

    if args.n_profiles_check > 0:
        make_sure_the_profile_is_correct(n_motifs_list, n_seqs_list,
                                         decompressed_profiles_array,
                                         args.are_seeds_degenerate,
                                         n_to_check = args.n_profiles_check,
                                         do_print=True)

    fold_instances()



if __name__ == "__main__":
    main()
