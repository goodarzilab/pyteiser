import numpy as np
import glob_var

class s_motif:

    def __init__(self, length):
        self.sequences = np.ones(shape=length, dtype=np.uint8)
        self.structures = np.ones(shape=length, dtype=np.uint8)
        self.num_phrases = length
        self.linear_length = length

    def print_sequence(self):
        string_to_print = ''.join([glob_var._char_to_nt_mapping[x] for x in self.sequences])
        print(string_to_print)

    def create_structure(self, stem_length, loop_length):
        # creates a stem-loop structure with the specified stem length and loop length
        # all bases will be U
        self.num_phrases = stem_length + loop_length
        self.linear_length = (stem_length * 2) + loop_length



# Pseudocode of lcl_create_motifs:
# the list we are putting all the motifs to is named motifs
# num_phrases = stem_length + loop_length
# lcl_initialize_sequence: initialize an array of NUCBITs of the length num_phrases and fill it up with 1s (U)
# Loop begins
# lcl_count_informative_bases: counts the number of bases that are not N and put it into num_inf_bases
# exit point: if( num_inf_bases < min_inf_bases or num_inf_bases < max_inf_bases ): skip this motif
# create_structure
# lcl_copy_sequence
# I = calcluate_information(motif)
# exit point: if( I < minI || I > maxI ): skip this motif
# if( num_motifs >= num_motifs_per_file ): flush buffer to file
# lcl_next_sequence( sequence, num_phrases )
# Loop ends
#
# Pseudocode of create_and_write_motifs:
# Loop 1 begins:
# for stem_length in range( min_stem_length, max_stem_length)
# Loop 2 begins:
# for loop_length in range( min_loop_length, max_loop_length)
# Loop 1 ends
# Loop 2 ends
#
#
#
#
#
#
#

#define _U				0x01 -> 01
#define _C				0x02 -> 010
#define _G				0x04 -> 0100
#define _A				0x08 -> 01000
#define _N				0x0F -> 01111
#define NUM_LETTERS		4
#define LAST_LETTER		_A


