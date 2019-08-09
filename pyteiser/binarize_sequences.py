import argparse

import IO


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seedfile", type=str)
    parser.add_argument("--rna_fastafile", type=str)
    parser.add_argument("--outfile", type=str)


    parser.set_defaults(
        rna_fastafile = '/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/narrow_down_transcripts_list/Gencode_v28_GTEx_expressed_transcripts_fasta/utr_3_fasta/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.txt',
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
    )

    args = parser.parse_args()

    return args



def main():
    args = handler()
    sequences_dict, seqs_order = IO.read_fasta(args.rna_fastafile)
    print("%d motifs have been loaded" % len(motifs_list))
    print(len(sequences_dict), len(seqs_order))




if __name__ == "__main__":
    main()