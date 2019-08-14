import argparse
import numpy as np
import timeit

import os
import sys

# to make sure relative import works
# for a detailed explanation, see test_matchmaker.py

current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)

import pyteiser.glob_var as glob_var
import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.matchmaker as matchmaker



def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seedfile", type=str)
    parser.add_argument("--rna_fastafile", type=str)
    parser.add_argument("--outfile", type=str)


    parser.set_defaults(
        rna_fastafile = '/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/narrow_down_transcripts_list/Gencode_v28_GTEx_expressed_transcripts_fasta/utr_3_fasta/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.txt',
        # rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        rna_bin_file='/Users/student/Documents/hani/temp/temp_bins/test_bin.bin',
    )

    args = parser.parse_args()

    return args


def main():
    args = handler()
    with open(args.rna_bin_file, 'rb') as rb:
        bitstring = rb.read()
        seq_objects_dict, seq_objects_order = IO.decompress_named_sequences(bitstring)
        full_string = IO.write_named_seq_to_fasta(seq_objects_dict, seq_objects_order)
    with open(args.rna_fastafile, 'r') as rf:
        full_fasta_string = rf.read()

    with open("/Users/student/Documents/hani/temp/temp_fasta/1.txt",'w') as wf:
        wf.write(full_string)

    full_fasta_string_Us = full_fasta_string.replace('T','U').replace('ENSU','ENST')
    assert(len(full_string) == len(full_fasta_string_Us))
    #
    # print(full_string[0:200])
    # print(full_fasta_string_Us[0:200])
    assert(full_string == full_fasta_string_Us)


if __name__ == "__main__":
    main()
