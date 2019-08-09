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

        current_motif = structures.s_motif(stem_length, loop_length)
        current_motif.sequence = curr_sequence
        current_motif.structure = curr_structure
        current_motif.compress()

        # current_motif.print_sequence()
        # current_motif.print_structure()

        assert(md5_checksum == current_motif.md5)

        motifs_list.append(current_motif)
        current_spot += 2 + 2*full_length + 16

    return motifs_list


def read_motif_file(inp_file):
    with open(inp_file, 'rb') as rf:
        full_bitstring = rf.read()
        motifs_list = decompress_motifs_from_bitstring(full_bitstring)
    return motifs_list


def read_fasta(infile):
    tr_dict_loc = {}
    seqs_order = []
    with open(infile, 'r') as f:
        split_string = f.read().split('>')
        for entry in split_string:
            if entry == '':
                continue
            seq_start = entry.find('\n')
            annotation = entry[:seq_start]
            sequence_string = entry[seq_start + 1:].replace('\n', '')
            current_sequence = structures.s_sequence(len(sequence_string))
            current_sequence.from_sequence(sequence_string)
            tr_dict_loc[annotation] = current_sequence
            seqs_order.append(annotation)

    return tr_dict_loc, seqs_order


def compress_named_sequences(seq_objects_dict, seqs_order):
    seq_batch_byte_string = b''

    for name in seqs_order:
        full_length_name = name.ljust(glob_var.MAX_SEQ_NAME_LENGTH)
        name_in_bytes = full_length_name.encode('utf-8')
        assert(len(name_in_bytes) == glob_var.MAX_SEQ_NAME_LENGTH)
        seq_batch_byte_string += name_in_bytes

        seq_objects_dict[name].compress()
        seq_batch_byte_string += seq_objects_dict[name].bytestring

    return seq_batch_byte_string
