import numpy as np
import pandas as pd
import argparse

import IO
import matchmaker
import type_conversions

def handler(raw_args = None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna_fastafile", help="fasta file with RNA sequences", type=str)
    parser.add_argument("--rna_bin_file", help="", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)

    parser.add_argument("--combined_seeds_filename", help="", type=str)
    parser.add_argument("--combined_profiles_filename", help="", type=str)
    parser.add_argument("--combined_MI_pv_zscores_filename", help="", type=str)
    parser.add_argument("--combined_robustness_filename", help="", type=str)

    parser.add_argument("--out_info_table", help="output file with statistics for the discovered seeds", type=str)
    parser.add_argument("--out_matches_table", help="output file with sequences of the matches for the discovered seeds", type=str)

    parser.add_argument("--indices_mode", help="compression in the index mode", type=bool)

    parser.set_defaults(
        rna_fastafile='/Users/student/Documents/hani/programs/pyteiser/data/tutorial_example_files/test_seqs.fa',
        rna_bin_file='/Users/student/Documents/hani/programs/pyteiser/data/temp/inputs/rna.bin',
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/temp/inputs/exp_mask.bin',
        combined_seeds_filename='/Users/student/Documents/hani/programs/pyteiser/data/temp/interm/optimized_seeds_1.bin',
        combined_profiles_filename='/Users/student/Documents/hani/programs/pyteiser/data/temp/interm/optimized_profiles_1.bin',
        combined_MI_pv_zscores_filename='/Users/student/Documents/hani/programs/pyteiser/data/temp/interm/optimized_MI_pv_zscores_1.bin',
        combined_robustness_filename='/Users/student/Documents/hani/programs/pyteiser/data/temp/interm/robustness_array_1.bin',
        out_info_table='/Users/student/Documents/hani/programs/pyteiser/data/temp/out/info.csv',
        out_matches_table='/Users/student/Documents/hani/programs/pyteiser/data/temp/out/matches.csv',

        indices_mode=True,

        # combined_seeds_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/tarbp2/seed_optimized_100k_tarbp2_utrs_10k.bin',
        # combined_profiles_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/tarbp2/profiles_optimized_100k_tarbp2_utrs_10k.bin',
        # combined_MI_pv_zscores_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/tarbp2/seed_characteristics_optimized_100k_tarbp2_utrs_10k.bin',
        # combined_robustness_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/tarbp2/robustness_optimized_100k_tarbp2_utrs_10k.bin',
        #
        # indices_mode=False,

        # combined_seeds_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/snrnpa1/seed_optimized_100k_snrnpa1_10k.bin',
        # combined_profiles_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/snrnpa1/profiles_optimized_100k_snrnpa1_10k.bin',
        # combined_MI_pv_zscores_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/snrnpa1/seed_characteristics_optimized_100k_snrnpa1_10k.bin',
        # combined_robustness_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/snrnpa1/robustness_optimized_100k_snrnpa1_10k.bin',
    )

    args = parser.parse_args(raw_args)

    return args


def read_seeds_and_characteristics(args):
    seeds_optimized = IO.read_motif_file(args.combined_seeds_filename)
    profiles_optimized = IO.unpack_profiles_file(args.combined_profiles_filename, args.indices_mode)
    seed_charact_array = IO.read_np_array(args.combined_MI_pv_zscores_filename, np.dtype('float64'))
    robustness_array = IO.read_np_array(args.combined_robustness_filename, np.dtype('bool'))

    return seeds_optimized, profiles_optimized, \
           seed_charact_array, robustness_array


def get_n_seqs_list(rna_bin_filename):
    seqs_dict, seqs_order = IO.read_rna_bin_file(rna_bin_filename)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)
    return n_seqs_list


def make_sequences_dataframe(fasta_filename, exp_mask_file):
    tr_dict_loc, seqs_order = IO.read_fasta_no_compression(fasta_filename)
    index_array, values_array = IO.unpack_mask_file(exp_mask_file)
    present_seqs_order = [x for i, x in enumerate(seqs_order) if index_array[i]]
    sub_dict_seq = {x : tr_dict_loc[x] for x in present_seqs_order}
    seq_df = pd.DataFrame.from_dict(sub_dict_seq, orient = 'index')
    seq_df = seq_df.rename({0 : 'sequence'}, axis = 1)
    return seq_df


def generate_report_csv(seeds_optimized, profiles_optimized,
                          seed_charact_array, robustness_array):
    MI_values_array = seed_charact_array[:, 0]
    seed_indices_sorted = np.argsort(MI_values_array)[::-1]

    report_csv = pd.DataFrame(columns = ['seed_id', 'sequence', 'structure',
                                         'MI', 'pvalue', 'zscore', 'robust', 'sequences_bound'],
                              index = np.arange(seed_indices_sorted.shape[0]))

    for index, i in enumerate(seed_indices_sorted):
        report_csv.at[index, 'seed_id'] = "%d" % i
        report_csv.at[index, 'sequence'] = seeds_optimized[i].print_linear_sequence(return_string = True)
        report_csv.at[index, 'structure'] = seeds_optimized[i].print_linear_structure(return_string = True)
        report_csv.at[index, 'MI'] = seed_charact_array[i,0]
        report_csv.at[index, 'pvalue'] = seed_charact_array[i,1]
        report_csv.at[index, 'zscore'] = seed_charact_array[i,2]
        report_csv.at[index, 'robust'] = robustness_array[i]
        report_csv.at[index, 'sequences_bound'] = profiles_optimized[i].sum()

    return report_csv


def add_matches_columns(seeds_optimized,
                        seed_charact_array,
                        n_seqs_list,
                        inp_df,
                        seq_column_name = 'sequence'):
    MI_values_array = seed_charact_array[:, 0]
    seed_indices_sorted = np.argsort(MI_values_array)[::-1]

    out_df = inp_df.copy()

    for index, i in enumerate(seed_indices_sorted):
        current_seed = seeds_optimized[i]
        current_seed.print_linear_sequence()
        current_seed.print_linear_structure()
        column_name = "matches_seed_%d" % (i)
        curr_matches_list = []
        for index in range(inp_df.shape[0]):
            sequence = inp_df.iloc[index][seq_column_name]
            all_instances = matchmaker.find_all_motif_instances(
                                                     type_conversions.w_to_n_motif(current_seed),
                                                     n_seqs_list[index],
                                                     is_degenerate = True)
            matches = []
            for inst in all_instances:
                curr_match = sequence[inst : inst + current_seed.linear_length]
                matches.append(curr_match)
            matches_string = "; ".join(matches)
            curr_matches_list.append(matches_string)
        out_df[column_name] = curr_matches_list
    return out_df




def main(raw_args = None):
    args = handler(raw_args)

    seeds_optimized, profiles_optimized, \
    seed_charact_array, robustness_array = read_seeds_and_characteristics(args)
    n_seqs_list = get_n_seqs_list(args.rna_bin_file)


    report_csv = generate_report_csv(seeds_optimized, profiles_optimized,
                                     seed_charact_array, robustness_array)

    seq_df = make_sequences_dataframe(args.rna_fastafile, args.exp_mask_file)
    matches_df = add_matches_columns(seeds_optimized,
                        seed_charact_array,
                        n_seqs_list,
                        seq_df)

    report_csv.to_csv(args.out_info_table, index = False, sep = '\t')
    matches_df.to_csv(args.out_matches_table, sep='\t')


if __name__ == "__main__":
    main()