import os
import sys
import numpy as np
import argparse
import math
import subprocess


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
        temp_folder = '/Users/student/Documents/hani/programs/pyteiser/data/temp_rnafold',
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
    global glob_var

    import structures
    import IO
    import type_conversions
    import glob_var


def calculate_window_length(args, do_print = True):
    window_length = args.max_stem_length * 2 + args.max_loop_length + \
        args.elongation_buffer_length + args.surrounding_buffer_length * 2
    if do_print:
       print("Folding sequences from the reference transcriptome "
             "with a rolling window of %d " % (window_length))
    return window_length


def generate_random_file_name(length_name = 7):
    rand_int_array = np.random.randint(10, size=length_name, dtype=np.uint8)
    rand_str_array = [str(x) for x in rand_int_array]
    rand_str = "".join(rand_str_array)
    return rand_str


def encode_folded_structure(dot_string, sequence_length):
    encoded_structure = np.zeros(sequence_length, dtype=np.uint8)

    for i in range(sequence_length):
        current_symbol = dot_string[i]
        current_encoding = glob_var._char_to_extended_structure[current_symbol]
        encoded_structure[i] = current_encoding

    return encoded_structure


def parse_RNAfold_output(infile):
    with open(infile, 'r') as rf:
        full_string = rf.read()
        splitted_string = full_string.split('\n')
        sequence_itself = splitted_string[0]
        sequence_length = len(sequence_itself)
        folded_structure = splitted_string[1][0 : sequence_length]
        encoded_structure = encode_folded_structure(folded_structure, sequence_length)

    return encoded_structure


def call_RNAfold(curr_sequence, args):
    rand_str = generate_random_file_name()
    fasta_sequence_filename = os.path.join(args.temp_folder, rand_str + '.fa')
    RNAfold_output_filename = os.path.join(args.temp_folder, rand_str + '_folded.fa')

    with open(fasta_sequence_filename, 'w') as wf:
        wf.write(curr_sequence)

    command_to_run = "RNAfold --noPS --noLP < %s" % (fasta_sequence_filename)
    with open(os.path.join(RNAfold_output_filename), 'w') as fout:
        subprocess.call(command_to_run, stdout=fout, shell=True)

    encoded_structure = parse_RNAfold_output(RNAfold_output_filename)

    os.remove(fasta_sequence_filename)
    os.remove(RNAfold_output_filename)

    return encoded_structure


def chunk_up_one_sequence(w_sequence, window_length, args):
    sequence_string = w_sequence.print(return_string = True)
    number_starting_points = len(sequence_string) - window_length + 1
    folded_structures_array = np.zeros((number_starting_points, window_length))

    for i in range(number_starting_points):
        current_sequence = sequence_string[i : i + window_length]
        encoded_structure = call_RNAfold(current_sequence, args)
        folded_structures_array[i] = encoded_structure
        #print(current_sequence)

        if i > 2:
            break

    print(folded_structures_array.shape)
    print(folded_structures_array)

    return folded_structures_array





def fold_all_sequences_wrapper(seqs_dict, seqs_order, window_length, args):
    for seq_name in seqs_order:
        folded_structures_array = chunk_up_one_sequence(seqs_dict[seq_name], window_length, args)
        current_bytestring = IO.write_np_array(folded_structures_array, out_filename = '', return_bytestring=False)

        break


def main():
    import_modules()
    args = handler()

    seqs_dict, seqs_order = IO.read_rna_bin_file(args.rna_bin_filename)
    window_length = calculate_window_length(args, do_print = True)
    fold_all_sequences_wrapper(seqs_dict, seqs_order, window_length, args)



if __name__ == "__main__":
    main()
