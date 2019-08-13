import numpy as np
import numba
import time

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

# I have tried really hard to improve performance of this step with numba
# the main problem is that I have a list of n_sequence objects and their size can vary
# therefore, I can't pass them to function as a numpy array with any of the standard formats
# I can mane a numpy array with an object dtyo (like dtype=structures.n_sequence) but Numba does not support it
# for more detailed explanation, see https://stackoverflow.com/questions/14639496/how-to-create-a-numpy-array-of-arbitrary-length-strings
# numba will deprecate standard python lists too
# there is also numba typed list structure (from numba.typed import List) but it's an experimental feature so far so I don't want to rely on it
# see here https://numba.pydata.org/numba-doc/dev/reference/pysupported.html
# so there is no way to pass a bunch of variable-sized sequence objects to numba in the way that would make the iterations faster


def calculate_profile_one_motif(motif, n_seqs_list):
    start_time = time.time()

    current_profile = structures.w_profile(len(n_seqs_list))
    for i, seq in enumerate(n_seqs_list):
        match = is_there_motif_instance(motif, seq)
        if match:
            current_profile.values[i] = True
    end_time = time.time()
    time_spent = end_time - start_time

    return current_profile, time_spent


def calculate_profiles_list_motifs(n_motifs_list, n_seqs_list,
                                   do_print=False):
    profiles_list = [0] * len(n_motifs_list)

    for i, motif in enumerate(n_motifs_list):
        current_profile, time_spent = calculate_profile_one_motif(motif, n_seqs_list)
        profiles_list[i] = current_profile.values

        if do_print:
            print("Motif number %d binds %d sequences. It took %.2f seconds"
                  % (i, current_profile.sum(), time_spent))

    return profiles_list
