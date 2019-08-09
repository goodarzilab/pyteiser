import numpy as np

import glob_var
import structures
import IO



def match_motif_seq(motif, sequence, ind):
    left_index = ind
    right_index = left_index + motif.linear_length - 1

    for i in range(motif.length):
        if motif.structure[i] == glob_var._stem:
            if not sequence.nt_is_a(left_index, motif.sequence[i]) or not sequence.is_paired(left_index, right_index):
                # either the sequence does not match or the left and right bases cannot form a Watson-Crick base pair
                return False
            left_index += 1
            right_index -= 1
        else:
            if not sequence.nt_is_a(left_index, motif.sequence[i]): # the sequence does not match
                return False
            left_index += 1
    return True




def is_there_motif_instance(motif, sequence):
    for i in range(sequence.length - motif.linear_length + 1):
        sequence_string = sequence.print(return_string = True)
        print(sequence_string[i : i + motif.linear_length])

        if match_motif_seq(motif, sequence, i):
            return True
    return False



def find_all_motif_instances(motif, sequence):
    motif_instances = []
    for i in range(sequence.length - motif.linear_length + 1):
        sequence_string = sequence.print(return_string = True)
        print(sequence_string[i : i + motif.linear_length])

        if match_motif_seq(motif, sequence, i):
            motif_instances.append(i)
    return motif_instances











def testing():
    test_seeds_file = '/Users/student/Documents/hani/temp/seeds_temp/python_generated_seeds/seeds_4-7_4-9_4-6_14-20_100_181.bin'
    test_seeds_list = IO.read_motif_file(test_seeds_file)
    seed_one = test_seeds_list[0]
    seed_one.print_sequence()
    seed_one.print_structure()

    test_sequence_string = 'UUUUUUUGACAACAAUUTGTCUUUUU'

    list_test_strings = ['UUUUUUUGACAACAAUUTGTCUUUUU',
                         'UUUUUUUGACAACAAUUTGTCUUUUU',
                         'UUUUUUUGACAACAAUUTGTCUUUUU']


    print("Searching for motif")
    for str in list_test_strings:
        test_sequence = structures.s_sequence(len(test_sequence_string))
        test_sequence.from_sequence(test_sequence_string)
        if is_there_motif_instance(seed_one, test_sequence):
            print("The sequence %s does contain the motif" % (str))







if __name__ == "__main__":
    testing()