import numpy as np
import argparse
import os

import glob_var
import structures
import IO
import matchmaker


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seedfile", type=str)
    parser.add_argument("--rna_fastafile", type=str)
    parser.add_argument("--outfile", type=str)


    parser.set_defaults(
        seedfile = '/Users/student/Documents/hani/temp/seeds_temp/python_generated_seeds/seeds_4-7_4-9_4-6_14-20_100_181.bin',
        rna_fastafile = '/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/narrow_down_transcripts_list/Gencode_v28_GTEx_expressed_transcripts_fasta/utr_3_fasta/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.txt',
        outfile = '/Users/student/Documents/hani/temp/profiles_temp/python_generated_profiles/profiles_4-7_4-9_4-6_14-20_100_181.bin'
    )

    args = parser.parse_args()

    return args



def main():
    args = handler()
    motifs_list = IO.read_motif_file(args.seedfile)
    sequences_dict, seqs_order = IO.read_fasta(args.rna_fastafile)
    print("%d motifs have been loaded" % len(motifs_list))
    print(len(sequences_dict), len(seqs_order))




if __name__ == "__main__":
    main()