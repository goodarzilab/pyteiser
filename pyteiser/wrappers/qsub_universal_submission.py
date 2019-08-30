import argparse
import random
import subprocess
import os
import sys


# The overall submission scheme:
# This (qsub_universal_submission.py) is the master submission script
# You call it with the arguments, specifying:
# 1) the Python script you want to submit
# 2) list of input files you want to run it on (more specifically, input/output folders, prefixes and list of file numbers)
# 3) list of submission arguments, like time and memory constraints
# 4) list of arguments that need to be passed to Python script you want to run
#
# This script then creates a qsub command of a standard form of
# qsub script_to_call [qsub options] [array of jobs t] [script parameters]
# and then submits it.
# Each individual script that is called through qsub then has to:
# 1) move to the working directory of the original script and one level up to import the rest of pyteiser
# 2) read the environmental variables to find out what task number was assigned to it
# 3) read the task number-to-input file mapping file to find out which file should it work on
# 4) read its own parameters
# 5) run
#
# Theoretically, I could submit each individual script with all the necessary qsub parameters avoiding such
# complicated scheme with an external wrapper
# The reason for doing things is that stupid SGE doesn't allow listing the indices of array jobs with comma
# (even though the documentation says it does). Therefore, if I have, say, files number 13,27,57 and I want to run
# my script on them, I have to either submit 3 jobs (bad) or create a separate file with mapping
# To avoid creating mapping files for each script individually, I use a system where a wrapper script does
# it automatically given the list of files that need to be processed


def handler():
    parser = argparse.ArgumentParser()

    #script to submit itself
    parser.add_argument("--script_to_sumbit", help="", type=str)

    # input-output: universal parameter names for all the scripts
    parser.add_argument("--input_folder", help="", type=str)
    parser.add_argument("--output_folder", help="", type=str)
    parser.add_argument("--input_file_prefix", help="", type=str)
    parser.add_argument("--output_file_prefix", help="", type=str)

    # list of input file indices (to create mapping to task numbers)
    parser.add_argument("--input_indices_list_file", help="input: list of indices of files to process", type=str)
    parser.add_argument("--mapping_task_ids_file", help="output: mapping of file indices to task ids", type=str)
    parser.add_argument("--", help="", type=str)
    parser.add_argument("--", help="", type=str)
    parser.add_argument("--", help="", type=str)

    # qsub submission parameters
    parser.add_argument("--python_binary", help="S parameter of qsub", type=str)
    parser.add_argument("--time_required", help="-l h_rt parameter of qsub", type=str)
    parser.add_argument("--mem_free", help="-l mem_free parameter of qsub", type=str)
    parser.add_argument("--mem_scratch", help="-l scratch parameter of qsub", type=str)
    parser.add_argument("--stderr_file", help="-e parameter of qsub", type=str)
    parser.add_argument("--stdout_file", help="-o parameter of qsub", type=str)
    parser.add_argument("--queue", help="-q parameter of qsub", type=str)
    parser.add_argument("--restart", help="-r parameter of qsub", type=str)

# script-specific arguments go into the -ac argument

    parser.set_defaults(
        incomplete_files_summary='/Users/student/Documents/hani/programs/pyteiser/data/testing_data/are_profiles_complete/profiles_tarbp2.txt',
        python_binary='/wynton/home/goodarzi/khorms/miniconda3/bin/python',
        time_required='60:00:00',
        mem_free='1G',
        mem_scratch='1G',
        stderr_file='/wynton/scratch/khorms/logs/test_stderr.txt',
        stdout_file='/wynton/scratch/khorms/logs/test_stdout.txt',
        queue='long.q',
        restart='yes',
    )

    args, unknown = parser.parse_known_args()

    return args, unknown


def create_mapping(args):
    with open(args.input_indices_list_file, 'r') as rf:
        full_string = rf.read().rstrip()
    indices_list_str = full_string.split(', ')
    indices_list = [int(x) for x in indices_list_str]

    random_int_list = [str(random.randint(0, 9)) for p in range(7)]
    random_int_string = "".join(random_int_list)
    unique_masking_file = "%s_%s" % (args.mapping_task_ids_file, random_int_string)
    with open(unique_masking_file, 'w') as wf:
        for task_id, index in enumerate(indices_list):
            current_string = "%d\t%d\n" % (task_id, index)
            wf.write(current_string)



def unknown_args_to_string(unknown_args):
    return " ".join(unknown_args)



def main():
    args, unknown_args = handler()






if __name__ == "__main__":
    main()