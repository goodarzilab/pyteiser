import numpy as np
import argparse
import hashlib

import glob_var
import structures



def decompress_motifs_from_bitstring(bitstring):
    motifs_list = []
    total_length = len(bitstring)
    current_spot = 0

    while current_spot < total_length:
        stem_length = bitstring[current_spot]
        loop_length = bitstring[current_spot + 1]
        full_length = stem_length + loop_length


        curr_sequence = np.frombuffer(bitstring[current_spot + 2 : current_spot + 2 + full_length], dtype=np.uint8)
        curr_structure = np.frombuffer(bitstring[current_spot + 2 + full_length :
                                        current_spot + 2 + 2*full_length], dtype=np.uint8)
        md5_checksum = bitstring[current_spot + 2 + 2*full_length: current_spot + 2 + 2*full_length + 16]

        current_motif = structures.w_motif(stem_length, loop_length)
        current_motif.sequence = curr_sequence
        current_motif.structure = curr_structure
        current_motif.compress()

        # current_motif.print_sequence()
        # current_motif.print_structure()

        assert(md5_checksum == current_motif.md5)

        motifs_list.append(current_motif)
        current_spot += 2 + 2*full_length + 16

    return motifs_list

def read_rna_bin_file(inp_file):
    with open(inp_file, 'rb') as rb:
        bitstring = rb.read()
        seq_objects_dict, seq_objects_order = decompress_named_sequences(bitstring)
    return seq_objects_dict, seq_objects_order


def read_motif_file(inp_file):
    with open(inp_file, 'rb') as rf:
        full_bitstring = rf.read()
        motifs_list = decompress_motifs_from_bitstring(full_bitstring)
    return motifs_list


def read_fasta(infile, do_print = False, how_often_print = 1000):
    tr_dict_loc = {}
    seqs_order = []
    with open(infile, 'r') as f:
        split_string = f.read().split('>')
        for ind, entry in enumerate(split_string):
            if entry == '':
                continue
            seq_start = entry.find('\n')
            annotation = entry[:seq_start]
            sequence_string = entry[seq_start + 1:].replace('\n', '')
            current_sequence = structures.w_sequence(len(sequence_string))
            current_sequence.from_sequence(sequence_string)
            tr_dict_loc[annotation] = current_sequence
            seqs_order.append(annotation)
            if ind % how_often_print == 0:
                if do_print:
                    print("Read sequence number ", ind)

    return tr_dict_loc, seqs_order


def compress_named_sequences(seq_objects_dict, seqs_order,
                             do_print=False, how_often_print=1000):

    seq_batch_byte_list = []

    for ind, name in enumerate(seqs_order):
        current_byte_string = b''
        full_length_name = name.ljust(glob_var.MAX_SEQ_NAME_LENGTH)
        name_in_bytes = full_length_name.encode('utf-8')
        assert(len(name_in_bytes) == glob_var.MAX_SEQ_NAME_LENGTH)
        current_byte_string += name_in_bytes

        seq_objects_dict[name].compress()

        current_byte_string += seq_objects_dict[name].bytestring
        seq_batch_byte_list.append(current_byte_string)
        if ind % how_often_print == 0:
            if do_print:
                print("Compressed sequence number ", ind)
    seq_batch_byte_string = b''.join(seq_batch_byte_list)

    return seq_batch_byte_string


def decompress_named_sequences(bitstring,
                               do_print=False, how_often_print=10000):
    seq_objects_dict = {}
    seq_objects_order = []

    total_length = len(bitstring)
    current_spot = 0
    counter = 0

    while current_spot < total_length:
        name_in_bytes = bitstring[current_spot : current_spot + glob_var.MAX_SEQ_NAME_LENGTH]
        length_bitstring = bitstring[current_spot + glob_var.MAX_SEQ_NAME_LENGTH :
                                        current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4]
        seq_length_np = np.frombuffer(length_bitstring, dtype=np.uint32)
        seq_length = seq_length_np[0]

        sequence_bitstring = bitstring[current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4 :
                                        current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4 + seq_length]
        md5_bitstring = bitstring[current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4 + seq_length :
                                        current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4 + seq_length + 16]

        current_spot += glob_var.MAX_SEQ_NAME_LENGTH + 4 + seq_length + 16

        full_length_name = name_in_bytes.decode('utf-8')
        name = full_length_name.rstrip()

        current_sequence = structures.w_sequence(seq_length)
        current_sequence.nts = np.frombuffer(sequence_bitstring, dtype = np.uint8)
        current_sequence.compress()


        assert (md5_bitstring == current_sequence.md5)

        seq_objects_dict[name] = current_sequence
        seq_objects_order.append(name)

        counter += 1
        if counter % how_often_print == 0:
            if do_print:
                print("Compressed sequence number ", counter)

    return seq_objects_dict, seq_objects_order


def write_named_seq_to_fasta(seq_objects_dict, seq_objects_order):
    strings_list = []
    for ind, name in enumerate(seq_objects_order):
        seq_string = seq_objects_dict[name].print_sequence(return_string = True)
        strings_list.append(">%s\n%s\n" % (name, seq_string))
    full_string = ''.join(strings_list)
    return full_string










