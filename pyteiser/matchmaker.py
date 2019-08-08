import numpy as np

import glob_var
import structures
import IO




def find_motif_instance():
    pass


def testing():
    test_seeds_file = '/Users/student/Documents/hani/temp/seeds_temp/python_generated_seeds/seeds_4-7_4-9_4-6_14-20_100_57.bin'
    test_seeds_list = IO.read_motif_file(test_seeds_file)
    seed_one = test_seeds_list[0]
    seed_one.print_sequence()
    seed_one.print_structure()
























if __name__ == "__main__":
    testing()