import argparse
import os
import sys

current_script_path = sys.argv[0]
sys.path.append("/wynton/home/goodarzi/khorms/pyteiser_root/pyteiser/pyteiser")

import glob_var


def get_potential_complementary_nts(k):
    parcitipants = glob_var._degenerate_nts_mapping[k]
    set_reversed = set()
    for pr in parcitipants:
        reversed_nts = glob_var._complementary_nt_sets_dict[pr]
        set_reversed = set_reversed.union(reversed_nts)

    return set_reversed


def find_corresponding_deg_nt(k, set_reversed, current_nt_char, compl_nt_dict, compl_nt_char_dict):
    for second_nt in glob_var._degenerate_nts_mapping:
        potential_set = glob_var._degenerate_nts_mapping[second_nt]
        if set_reversed == potential_set:
            compl_nt = second_nt
            compl_nt_char = glob_var._char_to_nt_mapping[compl_nt]

            compl_nt_dict[k] = compl_nt
            compl_nt_char_dict[current_nt_char] = compl_nt_char

    # in case if all nucleotides can be complementary to this one
    if len(set_reversed) == 4:
        compl_nt = glob_var._N
        compl_nt_char = glob_var._char_to_nt_mapping[compl_nt]

        compl_nt_dict[k] = compl_nt
        compl_nt_char_dict[current_nt_char] = compl_nt_char


    return compl_nt_dict, compl_nt_char_dict



def print_compl_dict(do_print = True):
    compl_nt_dict = {}
    compl_nt_char_dict = {}

    for k in glob_var.NT_LIST:
        current_nt_char = glob_var._char_to_nt_mapping[k]

        set_reversed = get_potential_complementary_nts(k)
        find_corresponding_deg_nt(k, set_reversed, current_nt_char, compl_nt_dict, compl_nt_char_dict)

    if do_print:
        print("Printing the complementary disctionary: ")
        for k in glob_var.NT_LIST:
            current_nt_char = glob_var._char_to_nt_mapping[k]
            compl_nt_char = glob_var._char_to_nt_mapping[compl_nt_dict[k]]

            print("_%s: _%s," % (current_nt_char, compl_nt_char))





if __name__ == "__main__":
    print_compl_dict()


