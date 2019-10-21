import numpy as np
import os
import sys

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import glob_var


def complementary_degenerate(deg_nt):
    deg_nt_char = glob_var._char_to_nt_mapping[deg_nt]
    participants = glob_var._degenerate_nts_mapping[deg_nt_char]
    print(participants)



def get_linear_sequence(w_motif):
    w_motif.print()

    full_sequence = np.zeros(w_motif.linear_length)

    for i in range(w_motif.length):
        left_index = i
        right_index = left_index + w_motif.linear_length - 1

        current_nt = w_motif.sequence[left_index]

        current_nt_char = glob_var._char_to_nt_mapping[current_nt]

        print(w_motif.structure[i], w_motif.sequence[i], current_nt_char)
        complementary_degenerate(current_nt)

        # full_sequence[left_index] = w_motif.sequence[left_index]
        # full_sequence[right_index] = complementary_degenerate(current_nt)





def generate_PWM_from_motif(inp_w_motif):
    pwm = np.zeros((inp_w_motif.linear_length, 4))

    labels = glob_var.PWM_LABELS

    get_linear_sequence(inp_w_motif)



    return pwm
# _degenerate_nts_mapping
