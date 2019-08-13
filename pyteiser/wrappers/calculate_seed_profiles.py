import numpy as np
import argparse
import sys
import os
import numba
import time

import glob_var
import structures
import IO
import matchmaker
import type_conversions


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seedfile", type=str)
    parser.add_argument("--rna_bin_file", type=str)
    parser.add_argument("--outfile", type=str)


    parser.set_defaults(
        seedfile = '/Users/student/Documents/hani/temp/seeds_temp/python_generated_seeds/seeds_4-7_4-9_4-6_14-20_100_181.bin',
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        outfile = '/Users/student/Documents/hani/temp/profiles_temp/python_generated_profiles/profiles_4-7_4-9_4-6_14-20_100_181.bin',
    )

    args = parser.parse_args()

    return args






def prepare_lists_for_calculations(args):
    seqs_dict, seqs_order = IO.read_rna_bin_file(args.rna_bin_file)
    w_motifs_list = IO.read_motif_file(args.seedfile)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]
    n_motifs_list = type_conversions.w_to_n_motifs_list(w_motifs_list)
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)
    matchmaker.calculate_profiles_list_motifs(n_motifs_list, n_seqs_list, do_print = True)


def main():
    args = handler()
    prepare_lists_for_calculations(args)


if __name__ == "__main__":
    main()