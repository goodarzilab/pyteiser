# run from /avicenna/khorms/projects/pyteiser_root/pyteiser/random_scripts/run_read_print_motifs.py

import os
import subprocess
import argparse


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--infolder", type=str)
    parser.add_argument("--outfolder", type=str)
    parser.add_argument("--prefix", type=str)
    parser.add_argument("--number_start", type=int) # which file number to start the transferring with
    parser.add_argument("--number_stop", type=int) # which file number to end the transferring with
    parser.add_argument("--c_script_path", type=str)

    parser.set_defaults(
        infolder='/avicenna/khorms/projects/pyteiser_root/debugging/test_generator/dump_c',
        outfolder = '/avicenna/khorms/projects/pyteiser_root/debugging/seeds_recoding/printed_text',
        prefix = 'seeds_c',
        number_start = 2214,
        number_stop = 2277,
        c_script_path='/avicenna/khorms/projects/iTEISER/c_programs_new/read_print_motifs',
    )

    args = parser.parse_args()

    return args



def generate_run_script(args):
    for i in range(args.number_start, args.number_stop + 1):
        c_seed_filename = os.path.join(args.infolder, "%s.%d.bin" % (args.prefix, i))
        text_seed_filename = os.path.join(args.outfolder, "%s.%d.txt" % (args.prefix, i))
        command_to_run = '%s -seedfile %s  -outfile %s' % (args.c_script_path, c_seed_filename, text_seed_filename)
        subprocess.call(command_to_run, shell=True)



def main():
    args = handler()
    generate_run_script(args)



if __name__ == "__main__":
    main()