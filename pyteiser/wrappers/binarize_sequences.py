import argparse
import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import IO


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna_fastafile", type=str)
    parser.add_argument("--rna_bin_file", type=str)


    parser.set_defaults(
        rna_fastafile = '/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/narrow_down_transcripts_list/Gencode_v28_GTEx_expressed_transcripts_fasta/utr_3_fasta/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.txt',
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        # rna_bin_file='/Users/student/Documents/hani/temp/temp_bins/test_bin.bin',
    )

    args = parser.parse_args()

    return args



def main():
    args = handler()
    sequences_dict, seqs_order = IO.read_fasta(args.rna_fastafile, do_print=True)
    seq_batch_byte_string = IO.compress_named_sequences(sequences_dict, seqs_order, do_print=True)
    with open(args.rna_bin_file, 'wb') as wf:
        wf.write(seq_batch_byte_string)



if __name__ == "__main__":
    main()