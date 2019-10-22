import numpy as np
import os
import sys

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import glob_var



def adjust_linear_length(w_motif):
    stem_count = np.count_nonzero(w_motif.structure == glob_var._stem)
    loop_count = np.count_nonzero(w_motif.structure == glob_var._loop)
    w_motif.linear_length = 2 * stem_count + loop_count



def get_linear_sequence(w_motif):
    w_motif.print()

    full_sequence = np.zeros(w_motif.linear_length)

    left_index = 0
    right_index = left_index + w_motif.linear_length - 1

    for i in range(w_motif.length):


        print(left_index, right_index)

        current_nt = w_motif.sequence[left_index]
        complementary_nt = glob_var._complementary_deg_nt_dict[current_nt]

        if w_motif.structure[i] == glob_var._stem:
            full_sequence[right_index] = complementary_nt
        full_sequence[left_index] = w_motif.sequence[left_index]

        left_index += 1
        right_index -= 1






def generate_PWM_from_motif(inp_w_motif):
    pwm = np.zeros((inp_w_motif.linear_length, 4))

    labels = glob_var.PWM_LABELS

    inp_w_motif.print()
    inp_w_motif.print_linear_sequence()

    # adjust_linear_length(inp_w_motif)
    # get_linear_sequence(inp_w_motif)



    return pwm
# _degenerate_nts_mapping
