import numpy as np
import argparse
import sys
import os
import numba

import glob_var
import structures
import IO
import matchmaker


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


def aaaa(seqs_dict, seqs_order, motifs_list, do_print = False,
         how_often_print=1000):
    for k1, motif in enumerate(motifs_list):
        current_sum = 0
        for k2, seq_name in enumerate(seqs_order):
            match = matchmaker.is_there_motif_instance(motif, seqs_dict[seq_name])
            if match:
                current_sum += 1
            if k2 % how_often_print == 0:
                if do_print:
                    print("Calculated matches for motif %d and sequence %d " % (k1, k2))
        curr_motif_string = motif.print_sequence(return_string = False)
        print("Motif %s binds %d sequences" % (curr_motif_string, current_sum))




def main():
    args = handler()
    seqs_dict, seqs_order = IO.read_rna_bin_file(args.rna_bin_file)
    motifs_list = IO.read_motif_file(args.seedfile)

    aaaa(seqs_dict, seqs_order, motifs_list, do_print = True)



if __name__ == "__main__":
    main()