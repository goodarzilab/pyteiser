import numpy as np
import glob_var
import structures
import argparse


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--min_stem_length", type=int)
    parser.add_argument("--max_stem_length", type=int)
    parser.add_argument("--min_loop_length", type=int)
    parser.add_argument("--max_loop_length", type=int)
    parser.add_argument("--print_sequences", type=str)
    parser.add_argument("--print_structures", type=str)
    # parser.add_argument("--", type=int)
    # parser.add_argument("--", type=int)
    # parser.add_argument("--", type=int)
    # parser.add_argument("--", type=int)


    parser.set_defaults(
        min_stem_length = 4,
        max_stem_length = 7,
        min_loop_length = 4,
        max_loop_length = 9,
        print_sequences = 'y',
        print_structures = 'n'
    )

    args = parser.parse_args()

    return args


def get_next_motif(motif, last_letter = glob_var._N):
    for i in range(motif.length):
        current_value = motif.sequences[i]
        if current_value < last_letter:
            motif.sequences[i] = current_value + glob_var._increment
            np.put(motif.sequences, ind=i, v=(current_value + glob_var._increment))
            # using numpy.put: https://docs.scipy.org/doc/numpy/reference/generated/numpy.put.html
            # alternatively, use numpy.left_shift: https://docs.scipy.org/doc/numpy/reference/generated/numpy.left_shift.html
            return True

        else:
            np.put(motif.sequences, ind=i, v=np.uint8(1))

    return False


# def testing():
#     toy_motif = structures.s_motif(length=8)
#     toy_motif.print_sequence()
#     for i in range(200):
#         get_next_motif(toy_motif)
#         toy_motif.print_sequence()


def create_motifs_fixed_length(stem_length, loop_length,
                               print_sequences, print_structures):
    toy_motif = structures.s_motif(stem_length, loop_length)

    if print_sequences == 'y':
        toy_motif.print_sequence()

    while get_next_motif(toy_motif):
        if print_sequences == 'y':
            toy_motif.print_sequence()



def create_all_motifs():
    args = handler()

    for stem_length in range(args.min_stem_length, args.max_stem_length + 1):
        for loop_length in range(args.min_loop_length, args.max_loop_length + 1):
            create_motifs_fixed_length(stem_length, loop_length,
                                       args.print_sequences, args.print_structures)



if __name__ == "__main__":
    create_all_motifs()
    # testing()



