import numpy as np
import argparse

import glob_var
import structures



def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--outfile", type=str)
    parser.add_argument("--min_stem_length", type=int)
    parser.add_argument("--max_stem_length", type=int)
    parser.add_argument("--min_loop_length", type=int)
    parser.add_argument("--max_loop_length", type=int)
    parser.add_argument("--print_sequences", type=str)
    parser.add_argument("--print_structures", type=str)
    parser.add_argument("--min_inf_bases", type=int)
    parser.add_argument("--max_inf_bases", type=int)
    parser.add_argument("--minI", type=int)
    parser.add_argument("--maxI", type=int)
    # parser.add_argument("--", type=int)
    # parser.add_argument("--", type=int)
    # parser.add_argument("--", type=int)
    # parser.add_argument("--", type=int)
    # parser.add_argument("--", type=int)
    # parser.add_argument("--", type=int)


    parser.set_defaults(
        outfile = '',
        min_stem_length = 4,
        max_stem_length = 7,
        min_loop_length = 4,
        max_loop_length = 9,
        min_inf_bases = 4,
        max_inf_bases = 6,
        minI = 14,
        maxI = 20,
        print_sequences = 'y',
        print_structures = 'y'
    )

    args = parser.parse_args()

    return args

def calculate_probability(motif):
  motif_probability = np.float64(0)

  for i in range(motif.stem_length):
      current_nt = motif.sequence[i]
      motif_probability += glob_var._paired_probabilities_dict[current_nt]

  for i in range(motif.stem_length, motif.length):
      current_nt = motif.sequence[i]
      if current_nt != glob_var._N:
          motif_probability += glob_var._paired_probabilities_dict['loop']

  motif_probability *= -1
  return motif_probability


def get_next_motif(motif, last_letter = glob_var._N):
    for i in range(motif.length):
        current_value = motif.sequence[i]
        if current_value < last_letter:
            motif.sequence[i] = current_value + glob_var._increment
            np.put(motif.sequence, ind=i, v=(current_value + glob_var._increment))
            # using numpy.put: https://docs.scipy.org/doc/numpy/reference/generated/numpy.put.html
            # alternatively, use numpy.left_shift: https://docs.scipy.org/doc/numpy/reference/generated/numpy.left_shift.html
            return True

        else:
            np.put(motif.sequence, ind=i, v=np.uint8(1))

    return False


# def testing():
#     toy_motif = structures.s_motif(length=8)
#     toy_motif.print_sequence()
#     for i in range(200):
#         get_next_motif(toy_motif)
#         toy_motif.print_sequence()

def count_informative_bases(motif):
    inf_bases = np.count_nonzero(glob_var._N != motif.sequence)
    return inf_bases


def check_motif_criteria(motif, args):
    # check for number of informative bases
    # and check for how much information does the motif carry

    num_inf_bases = count_informative_bases(motif)
    # print(num_inf_bases)
    if num_inf_bases < args.min_inf_bases or num_inf_bases > args.max_inf_bases:
        return False

    motif_probability = calculate_probability(motif)
    # print(motif_probability)
    if motif_probability < args.minI or motif_probability > args.maxI:
        return False

    return True


def create_motifs_fixed_length(stem_length, loop_length,
                               print_sequences, print_structures,
                               args):
    curr_motif = structures.s_motif(stem_length, loop_length)

    if print_sequences == 'y':
        curr_motif.print_sequence()
    if print_structures == 'y':
        curr_motif.print_structure()


    while get_next_motif(curr_motif):
        if not check_motif_criteria(curr_motif, args):
            continue


        if print_sequences == 'y':
            curr_motif.print_sequence()
        if print_structures == 'y':
            curr_motif.print_structure()





def create_all_motifs():
    args = handler()

    for stem_length in range(args.min_stem_length, args.max_stem_length + 1):
        for loop_length in range(args.min_loop_length, args.max_loop_length + 1):
            create_motifs_fixed_length(stem_length, loop_length,
                                       args.print_sequences, args.print_structures,
                                       args)



if __name__ == "__main__":
    create_all_motifs()



