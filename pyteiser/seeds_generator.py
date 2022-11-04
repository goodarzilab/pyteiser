import numpy as np
import argparse
import os
import sys

 try:
     from . import glob_var
 except:
     import glob_var

 try:
     from . import structures
 except:
     import structures



# this function generates all possible seeds with specified length of stem and loop
# to reduce the search space, it only keeps the seeds that have a pre-specified number of informative bases (non-Ns)
# and also that fall into a certain range of information content
# the ones that we keep are being compressed to a bit string (see structures module) with necessary MD5 checksum
# and then being written to series of files of specified size


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--outfolder", type=str)
    parser.add_argument("--prefix", type=str)
    parser.add_argument("--num_motifs_per_file", type=int)
    parser.add_argument("--min_stem_length", type=int)
    parser.add_argument("--max_stem_length", type=int)
    parser.add_argument("--min_loop_length", type=int)
    parser.add_argument("--max_loop_length", type=int)
    parser.add_argument("--print_sequences", type=str)
    parser.add_argument("--print_structures", type=str)
    parser.add_argument("--print_first_motif", type=str)
    parser.add_argument("--min_inf_bases", type=int)
    parser.add_argument("--max_inf_bases", type=int)
    parser.add_argument("--minI", type=int)
    parser.add_argument("--maxI", type=int)


    parser.set_defaults(
        outfolder = '/Users/student/Documents/hani/temp/seeds_temp/python_generated_seeds',
        prefix = 'seeds_4-7_4-9_4-6_14-20_100',
        num_motifs_per_file = 100,
        min_stem_length = 4,
        max_stem_length = 7,
        min_loop_length = 4,
        max_loop_length = 9,
        min_inf_bases = 4,
        max_inf_bases = 6,
        minI = 14,
        maxI = 20,
        print_sequences = 'n',
        print_structures = 'n',
        print_first_motif = 'n'
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


def create_motifs_fixed_length(stem_length, loop_length, motifs_counter, total_bitstring,
                               do_print_first_motif, do_print_sequences, do_print_structures,
                               args):
    curr_motif = structures.w_motif(stem_length, loop_length) # this is a motif with all U, we are skipping it anyway

    while get_next_motif(curr_motif):
        if not check_motif_criteria(curr_motif, args):
            continue
        motifs_counter += 1

        curr_motif.compress()
        total_bitstring += curr_motif.bytestring

        if do_print_first_motif:
            if (motifs_counter % args.num_motifs_per_file) == 1:
                print("In the file number %d the first motif is" % ((motifs_counter // args.num_motifs_per_file) + 1))
                curr_motif.print_sequence()
                sys.stdout.flush()

        if (motifs_counter % args.num_motifs_per_file) == 0:
            next_filename = os.path.join(args.outfolder, "%s_%d.bin" % (args.prefix, motifs_counter // args.num_motifs_per_file))
            with open(next_filename, 'wb') as wf:
                wf.write(total_bitstring)
            total_bitstring = b''

        if do_print_sequences:
            curr_motif.print_sequence()
        if do_print_structures:
            curr_motif.print_structure()

    return motifs_counter, total_bitstring



def process_printing_arguments(args):
    do_print_first_motif = False
    do_print_sequences = False
    do_print_structures = False

    if args.print_first_motif == 'y':
        do_print_first_motif = True
    if args.print_sequences == 'y':
        do_print_sequences = True
    if args.print_structures == 'y':
        do_print_structures = True

    return do_print_first_motif, do_print_sequences, do_print_structures



def motif_creator_outer_loop(args):
    do_print_first_motif, do_print_sequences, do_print_structures = process_printing_arguments(args)
    motifs_counter = 0
    total_bitstring = b''

    for stem_length in range(args.min_stem_length, args.max_stem_length + 1):
        for loop_length in range(args.min_loop_length, args.max_loop_length + 1):
            motifs_counter, total_bitstring = create_motifs_fixed_length(stem_length, loop_length, motifs_counter, total_bitstring,
                                                        do_print_first_motif, do_print_sequences, do_print_structures,
                                                        args)



def create_all_motifs():
    args = handler()
    motif_creator_outer_loop(args)



if __name__ == "__main__":
    create_all_motifs()



