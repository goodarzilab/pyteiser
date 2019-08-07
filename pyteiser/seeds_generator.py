import numpy as np
import glob_var
import structures

def get_next_motif(motif, last_letter = glob_var._N):
    for i in range(motif.linear_length):
        current_value = motif.sequences[i]
        if current_value < last_letter:
            motif.sequences[i] = current_value + glob_var._increment
            np.put(motif.sequences, ind=i, v=(current_value + glob_var._increment))
            # using numpy.put: https://docs.scipy.org/doc/numpy/reference/generated/numpy.put.html
            # alternatively, use numpy.left_shift: https://docs.scipy.org/doc/numpy/reference/generated/numpy.left_shift.html
            return

        else:
            np.put(motif.sequences, ind=i, v=np.uint8(1))


def testing():
    toy_motif = structures.s_motif(length=8)
    toy_motif.print_sequence()
    for i in range(200):
        get_next_motif(toy_motif)
        toy_motif.print_sequence()



if __name__ == "__main__":
    testing()


