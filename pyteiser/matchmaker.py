import numpy as np
import numba

import glob_var
import structures
import IO



@numba.jit(cache=True, nopython=True, nogil=True)
def match_motif_seq(n_motif, n_sequence, ind):
    # this function only works with n_motif and n_sequence classes,
    # not with w_motif and w_sequence

    left_index = ind
    right_index = left_index + n_motif.linear_length - 1

    for i in range(n_motif.length):
        if n_motif.structure[i] == glob_var._stem:
            if not n_sequence.nt_is_a(left_index, n_motif.sequence[i]) or not n_sequence.is_paired(left_index, right_index):
                # either the sequence does not match or the left and right bases cannot form a Watson-Crick base pair
                return False
            left_index += 1
            right_index -= 1
        else:
            if not n_sequence.nt_is_a(left_index, n_motif.sequence[i]): # the sequence does not match
                return False
            left_index += 1
    return True



@numba.jit(cache=True, nopython=True, nogil=True)
def is_there_motif_instance(n_motif, n_sequence):
    # this function only works with n_motif and n_sequence classes,
    # not with w_motif and w_sequence

    for i in range(n_sequence.length - n_motif.linear_length + 1):
        # sequence_string = sequence.print(return_string = True)
        # print(sequence_string[i : i + motif.linear_length])
        if match_motif_seq(n_motif, n_sequence, i):
            return True
    return False


@numba.jit(cache=True, nopython=True, nogil=True)
def find_all_motif_instances(n_motif, n_sequence):
    # this function only works with n_motif and n_sequence classes,
    # not with w_motif and w_sequence

    motif_instances = []
    for i in range(n_sequence.length - n_motif.linear_length + 1):
        # sequence_string = sequence.print(return_string = True)
        # print(sequence_string[i : i + motif.linear_length])
        if match_motif_seq(n_motif, n_sequence, i):
            motif_instances.append(i)
    return motif_instances



