import numpy as np
import argparse
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

    parser.add_argument("--incomplete_files_summary", help="list of incomplete profile files", type=str)
    parser.add_argument("--python_binary", help="S parameter of qsub", type=str)
    parser.add_argument("--time_required", help="-l h_rt parameter of qsub", type=str)
    parser.add_argument("--mem_free", help="-l mem_free parameter of qsub", type=str)
    parser.add_argument("--mem_scratch", help="-l scratch parameter of qsub", type=str)
    parser.add_argument("--stderr_file", help="-e parameter of qsub", type=str)
    parser.add_argument("--stdout_file", help="-o parameter of qsub", type=str)
    parser.add_argument("--queue", help="-q parameter of qsub", type=str)
    parser.add_argument("--restart", help="-r parameter of qsub", type=str)
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

    args = parser.parse_args()

    return args



def aaa():
    pass



def main():
    args = handler()






if __name__ == "__main__":
    main()