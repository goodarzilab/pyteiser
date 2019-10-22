import numpy as np
import pandas as pd
import seqlogo
import os
import sys

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import glob_var

import modify_seed
import type_conversions


def draw_pic():
    ACGT_frequencies = np.zeros((4, 4))
    ACGT_frequencies[0, 0] = 1
    ACGT_frequencies[1, 1] = 1
    ACGT_frequencies[1, 2] = 1
    ACGT_frequencies[2, 1] = 1
    ACGT_frequencies[2, 2] = 1
    ACGT_frequencies[2, 3] = 1
    ACGT_frequencies[3, 1] = 1
    ACGT_frequencies[3, 2] = 1
    ACGT_frequencies[3, 3] = 1
    ACGT_frequencies[3, 0] = 1

    ACGT_frequencies = ACGT_frequencies / ACGT_frequencies.sum(axis=1, keepdims=True)
    ACGT_frequencies_pd = pd.DataFrame(ACGT_frequencies)


    background_NA_dict = {nt: 0.25 for nt in 'ACGU'}
    background_NA_list = np.array(list(background_NA_dict.values()))

    ACGT_ppm_temp = seqlogo.Ppm(ACGT_frequencies_pd, alphabet_type='RNA', background=background_NA_list)
    ACGT_ppm_temp.idxmax = ACGT_frequencies_pd.idxmax
    ACGT_frequencies_ppm = seqlogo.Ppm(ACGT_ppm_temp, alphabet_type='RNA', alphabet='ACGU',
                                       background=background_NA_list)

    seqlogo.seqlogo(ACGT_frequencies_ppm, color_scheme='classic',
                    format='png', size='medium', filename='/Users/student/Desktop/a2.png')




def generate_PWM_from_motif(inp_w_motif, do_print = True):
    pwm = np.zeros((inp_w_motif.linear_length, 4))

    labels = glob_var.PWM_LABELS

    if do_print:
        print("Short motif representation:")
        inp_w_motif.print()
        print("Full motif representation:")
        inp_w_motif.print_linear()
        print()

    inp_n_motif = type_conversions.w_to_n_motif(inp_w_motif)
    m_seeds = modify_seed.elongate_motif(inp_n_motif)

    for seed in m_seeds:
        current_w_seed = type_conversions.n_to_w_motif(seed)
        current_w_seed.print()
        current_w_seed.print_linear()
        print()






    return pwm


draw_pic()