import numpy as np
import pandas as pd
import seqlogo
import os
import sys


current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

sys.path.append("/Users/student/Downloads/weblogo-master")
import weblogo as wl

import glob_var

import modify_seed
import type_conversions


def generate_frequencies(inp_w_motif):
    frequencies = np.zeros((inp_w_motif.linear_length, 4))
    inp_w_motif.get_linear_sequence()
    for i, nt in enumerate(inp_w_motif.linear_sequence):
        set_of_letters = glob_var._degenerate_nts_mapping[nt]
        for let in set_of_letters:
            column_index = glob_var._rna_alphabet_nt_mapping[let]
            frequencies[i, column_index] = 1
        if len(set_of_letters) == 0:
            frequencies[i, :] = np.array([1,1,1,1])

    return frequencies


def process_freq_for_seqlogo(frequencies):
    background_NA_dict = {nt: 0.25 for nt in glob_var._rna_alphabet_list}
    background_NA_list = np.array(list(background_NA_dict.values()))
    # create PFM object for seqlogo
    pfm = seqlogo.Pfm(frequencies, alphabet_type='RNA', background=background_NA_list)
    # compute missing data: PPM, PWM, information content etc
    # no need to do weird reformatting between PFM, PPM, PWM
    # see https://pypi.org/project/seqlogo/
    # don't look at "Generate some frequency data and convert to PWM" - it's misleading
    PM = seqlogo.CompletePm(pfm, alphabet=glob_var._rna_alphabet_string)

    return PM


def draw_weblogo(inp_w_motif, out_file):
    frequencies = generate_frequencies(inp_w_motif)
    plottable_matrix = process_freq_for_seqlogo(frequencies)

    seqlogo.seqlogo(plottable_matrix, color_scheme='classic',
                    format='eps', size='medium', filename=out_file)


# def print_motif(inp_w_motif, do_print = True):
#     if do_print:
#         print("Short motif representation:")
#         inp_w_motif.print()
#         print("Full motif representation:")
#         inp_w_motif.print_linear()
#         print()

def draw_secondary_structure(inp_w_motif):
    pass

