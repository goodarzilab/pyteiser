import numpy as np
import argparse
import hashlib

import glob_var
import structures



def decompress_motifs_from_bitstring(bitstring):
    total_length = len(bitstring)
    current_spot = 0

    while current_spot < total_length:
        stem_length = bitstring[current_spot]
        loop_length = bitstring[current_spot + 1]
        full_length = stem_length + loop_length

        print(stem_length, loop_length)

        curr_sequence = np.frombuffer(bitstring[current_spot + 2 : current_spot + 2 + full_length], dtype=np.uint8)
        curr_structure = np.frombuffer(bitstring[current_spot + 2 + full_length :
                                        current_spot + 2 + 2*full_length], dtype=np.uint8)
        md5_checksum = np.frombuffer(bitstring[current_spot + 2 + full_length:
                                        current_spot + 2 + 2*full_length + 16], dtype=np.uint8)

        current_motif = structures.s_motif(stem_length, loop_length)
        current_motif.sequence = curr_sequence
        current_motif.structure = curr_structure
        current_motif.compress()

        current_motif.print_sequence()
        current_motif.print_structure()

        assert(md5_checksum == current_motif.md5)



def read_motif_file(inp_file):
    with open(inp_file, 'rb') as rf:
        full_bitstring = rf.read()
        decompress_motifs_from_bitstring(full_bitstring)


read_motif_file('/Users/student/Documents/hani/temp/seeds_temp/python_generated_seeds/seeds_4-7_4-9_4-6_14-20_100_194.bin')