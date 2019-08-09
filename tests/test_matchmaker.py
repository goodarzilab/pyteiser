# see how it's implemented here: https://www.kennethreitz.org/essays/repository-structure-and-python

import os
import sys
import numpy as np
sys.path.insert(0, os.path.abspath('..'))

import pyteiser.glob_var as glob_var
import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.matchmaker as matchmaker




def test_nts():
    motif_1_string = "GNCANCNNUU"
    # motif_2 = "NNNCNACGUU"
    # motif_3 = "NGNGNGNCUU"

# 4 6 NNACGNNCUU
# 4 6 ANCCNNNUUU
# 4 6 NNCNANAGUU



def testing():
    test_seeds_file = '/Users/student/Documents/hani/temp/seeds_temp/python_generated_seeds/seeds_4-7_4-9_4-6_14-20_100_63.bin'
    test_seeds_list = IO.read_motif_file(test_seeds_file)
    seed_one = test_seeds_list[0]
    seed_two = test_seeds_list[1]
    seed_one.print_sequence()
    seed_one.print_structure()
    #
    # seed_one_cr = structures.s_motif(4,6)
    # seed_one_cr.from_string("GNCANCNNUU")
    # seed_one_cr.print_sequence()
    # seed_one_cr.print_structure()
    # if not seed_one.__eq__(seed_one_cr):
    #     print("not equal")
    # else:
    #     print("equal")
    #
    # if not seed_two.__eq__(seed_one_cr):
    #     print("not equal 2")
    # else:
    #     print("equal 2")
    #



        # test_sequence_string = 'UUUUUUUGACAACAAUUTGTCUUUUU'
    #
    # list_test_strings = ['UUUUUUUGACAACAAUUTGTCUUUUU',
    #                      'UUUUUUUGACAACAAUUTGTCUUUUU',
    #                      'UUUUUUUGACAACAAUUTGTCUUUUU']


    # print("Searching for motif")
    # for str in list_test_strings:
    #     test_sequence = structures.s_sequence(len(test_sequence_string))
    #     test_sequence.from_sequence(test_sequence_string)
    #     if matchmaker.is_there_motif_instance(seed_one, test_sequence):
    #         print("The sequence %s does contain the motif" % (str))







if __name__ == "__main__":
    testing()
