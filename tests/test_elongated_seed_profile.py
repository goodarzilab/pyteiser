import numpy as np
import argparse
import os
import sys
import math
import copy


current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)

import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.glob_var as glob_var
import pyteiser.matchmaker as matchmaker
import pyteiser.type_conversions as type_conversions
import pyteiser.MI as MI
import pyteiser.modify_seed as modify_seed


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rna_bin_file", help="referense transcriptome in binary format", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)

    parser.set_defaults(
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',
        nbins=15,
        number_example_matches_to_print = 5,
    )

    args = parser.parse_args()

    return args


def read_sequences(rna_bin_filename):
    seqs_dict, seqs_order = IO.read_rna_bin_file(rna_bin_filename)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)

    return n_seqs_list


def create_one_seed(do_print = True):
    test_motif_1 = structures.w_motif(5, 4)
    test_motif_1.from_string("BSHNVBCNU")
    test_motif_1.change_structure_position(0, glob_var._loop)
    if do_print:
        print("Test motif:")
        test_motif_1.print()
        test_motif_1.print_linear()
    n_test_motif = type_conversions.w_to_n_motif(test_motif_1)
    return n_test_motif


def test_elongated_seed(seqs_of_interest, discr_exp_profile, nbins, N, do_print=True):
    elong_seed = create_one_seed(do_print)
    current_profile, time_spent = matchmaker.calculate_profile_one_motif(elong_seed,
                                                                         seqs_of_interest,
                                                                         is_degenerate=True)
    matching_sequences = [seqs_of_interest[x] for x in range(current_profile.values.shape[0]) if current_profile.values[x]]
    curr_mi = MI.mut_info(current_profile.values, discr_exp_profile, x_bins=2, y_bins=nbins)
    if do_print:
        print(curr_mi)

    first_N_matching_sequences = matching_sequences[0:N]

    counter = 0

    for seq in first_N_matching_sequences:
        curr_matching_indices = matchmaker.find_all_motif_instances(elong_seed, seq, is_degenerate = True)
        for match_index in curr_matching_indices:
            counter += 1
            match_sequence = structures.w_sequence(elong_seed.linear_length)
            match_sequence.nts = seq.nts[match_index : match_index + elong_seed.linear_length]
            match_string = match_sequence.print(return_string = True)
            print("Match %d: %s" % (counter, match_string))




def main():
    args = handler()

    n_seqs_list = read_sequences(args.rna_bin_file)
    index_array, values_array = IO.unpack_mask_file(args.exp_mask_file)
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins = args.nbins)
    seqs_of_interest = [n_seqs_list[x] for x in range(index_array.shape[0]) if index_array[x]]

    test_elongated_seed(seqs_of_interest, discr_exp_profile, args.nbins, args.number_example_matches_to_print)







if __name__ == "__main__":
    main()


# test matches of a seed:
# BSHNVBCNU
# .<<<<....