# run from /avicenna/khorms/projects/pyteiser_root/pyteiser/random_scripts/text_motifs_to_pyteiser.py

import sys
import os
import subprocess
import argparse

pyteiser_home_path = "/avicenna/khorms/projects/pyteiser_root/pyteiser"
sys.path.append(pyteiser_home_path)

import pyteiser.structures as structures
import pyteiser.IO as IO


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--infolder", type=str)
    parser.add_argument("--outfolder", type=str)
    parser.add_argument("--prefix", type=str)
    parser.add_argument("--number_start", type=int) # which file number to start the transferring with
    parser.add_argument("--number_stop", type=int) # which file number to end the transferring with
    # parser.add_argument("--c_script_path", type=str)

    parser.set_defaults(
        infolder='/avicenna/khorms/projects/pyteiser_root/debugging/seeds_recoding/printed_text',
        outfolder = '/avicenna/khorms/projects/pyteiser_root/debugging/seeds_recoding/to_pytheiser',
        old_prefix = 'seeds_c.',
        new_prefix = 'seeds_4-7_4-9_4-6_14-20_30k_',
        number_start = 2214,
        number_stop = 2277,
        # c_script_path='/avicenna/khorms/projects/iTEISER/c_programs_new/read_print_motifs',
    )

    args = parser.parse_args()

    return args



def convert_motif_text_file(inp_filename, out_filename):
    current_sequence_str = ''
    current_structure_str = ''

    total_bitstring = b''

    with open(inp_filename, 'r') as rf:
        for i, line in enumerate(rf):
            stripped_line = line.rstrip()

            if i % 2 == 0:
                current_structure_str = stripped_line
            else:
                current_sequence_str = stripped_line

                stem_length = current_structure_str.index('.')
                loop_length = len(current_structure_str) - stem_length

                current_motif = structures.w_motif(stem_length, loop_length)
                current_motif.from_string(current_sequence_str)

                current_motif.compress()
                total_bitstring += current_motif.bytestring

    with open(out_filename, 'wb') as wf:
        wf.write(total_bitstring)

                




def generate_run_script(args):
    for i in range(args.number_start, args.number_stop + 1):
        text_seed_filename = os.path.join(args.infolder, "%s%d.txt" % (args.old_prefix, i))
        pyteiser_seed_filename = os.path.join(args.outfolder, "%s%d.bin" % (args.new_prefix, i))
        
        convert_motif_text_file(text_seed_filename, pyteiser_seed_filename)



def main():
    args = handler()
    generate_run_script(args)



if __name__ == "__main__":
    main()


# To test:

# pyteiser_home_path = "/Users/student/Documents/hani/programs/pyteiser"
# sys.path.append(pyteiser_home_path)

# import pyteiser.IO as IO

# or whatever other test file this script has generated
# test_seeds_filename = '/Users/student/Documents/hani/temp/converted_seeds/seeds_4-7_4-9_4-6_14-20_30k_2214.bin'

# w_motifs_list = IO.read_motif_file(test_seeds_filename)

# for i in range(8):
#     w_motifs_list[i].print()

