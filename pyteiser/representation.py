import numpy as np
import os
import sys

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import glob_var

import modify_seed
import type_conversions







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
