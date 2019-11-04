import os
import sys
import numpy as np
import argparse
import math


# this script folds sequences from the reference transcriptome
# most methods use dynamic programming algorithm that works in a cubic time
# it's fine when you work with short-range base pairing but it's very expensive to work with long base-paiting
# some other algorithms like LinearFold and co-transcriptional folding are designed to improve speed and quality
# of long base pairing prediction
# However, predictions of long-range basepairing don't work well anyway and we don't really care about them
# in the context of pyteiser anyway. Therefore, we are focusing on short range interactions for pyteiser
# so we fold each sequence in a rolling ~50-nt window and then compare the seed matches to the predicted folding
# we require the exact match of the two predicted structures


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna_bin_filename", type=str)
    parser.add_argument("--temp_folder", type=str)

    parser.add_argument("--max_stem_length", type=int)
    parser.add_argument("--max_loop_length", type=int)
    parser.add_argument("--elongation_buffer_length", type=int)
    parser.add_argument("--surrounding_buffer_length", type=int)



    parser.set_defaults(
        rna_bin_filename = '/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        max_stem_length = 7,
        max_loop_length = 9,
        elongation_buffer_length = 6,
        surrounding_buffer_length = 10,
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

    import structures
    import IO
    import type_conversions


def calculate_window_length(args, do_print = True):
    window_length = args.max_stem_length * 2 + args.max_loop_length + \
        args.elongation_buffer_length + args.surrounding_buffer_length * 2
    if do_print:
       print("Folding sequences from the reference transcriptome "
             "with a rolling window of %d " % (window_length))
    return window_length


def call_RNAfold(curr_sequence):
    # TODO: create random file name
    # TODO: run RNAfold
    fasta_sequence_file =
    command_to_run = "RNAfold --noPS --noLP < %s" % (full_name_confile, fasta_sequence_file)
    with open(os.path.join(current_outfile), 'w') as fout:
        subprocess.call(command_to_run, stdout=fout, shell=True)


def chunk_up_one_sequence(w_sequence, window_length, args):
    sequence_string = w_sequence.print(return_string = True)

    for i in range(len(sequence_string) - window_length + 1):
        current_sequence = sequence_string[i : i + window_length]
        print(current_sequence)




def fold_all_sequences_wrapper(seqs_dict, seqs_order, window_length, args):
    for seq_name in seqs_order:
        chunk_up_one_sequence(seqs_dict[seq_name], window_length, args)

        break


def main():
    import_modules()
    args = handler()

    seqs_dict, seqs_order = IO.read_rna_bin_file(args.rna_bin_filename)
    window_length = calculate_window_length(args, do_print = True)
    fold_all_sequences_wrapper(seqs_dict, seqs_order, window_length, args)



if __name__ == "__main__":
    main()
