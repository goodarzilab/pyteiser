#$ -S /wynton/home/goodarzi/khorms/miniconda3/bin/python
#$ -t 1-2277
#$ -l h_rt=90:00:00
#$ -l mem_free=100G
#$ -l scratch=100G
#$ -e /wynton/scratch/khorms/logs/test_stderr.txt
#$ -o /wynton/scratch/khorms/logs/test_stdout.txt
#$ -q long.q
#$ -r yes

# --------------------------------
# QSUB FROM ITS PARENT FOLDER ONLY
# SEE EXPLANATION BELOW
# --------------------------------

# for some reason qsub doesn't want to implement python scripts when I point to an inside-of-conda-environment python
# executable, unless it's the base environment.
# for example, if I want to run a python script through qsub in the pyteiser_env environment, I get the path to executable with
# import sys
# print(sys.executable)
# and it gives me /wynton/home/goodarzi/khorms/miniconda3/envs/pyteiser_env/bin/python
# however, if I put this executable in the first line of my script I get an error
# described here https://stackoverflow.com/questions/19292957/how-can-i-troubleshoot-python-could-not-find-platform-independent-libraries-pr
# so instead I have to point it to an executable in the base environment, namely /wynton/home/goodarzi/khorms/miniconda3/bin/python
# it seems like it has something to do with PYTHONPATH and PYTHONHOME environmental variables, but I don't know how to fix it in a qsub script
# see here https://stackoverflow.com/questions/19292957/how-can-i-troubleshoot-python-could-not-find-platform-independent-libraries-pr


import numpy as np
import argparse
import subprocess
import os
import sys


# when I submit an array job with qsub, SGE creates lots of individual scripts and runs them one by one
# these individual scripts run on the cluster nodes and they don't know where is the parent script stored
# the cluster nodes might not even have access to the submission directory. See here: https://stackoverflow.com/questions/37721931/sge-get-script-location-inside-the-script
# therefore, my hacky solution of python relative import problem doesn't work with qsub submissions, since
# it is based on the fact that master script knows where it has been submitted from.
# this script is built in the way where you have to run it from the same directory where it is located
# so usage is: (1) cd path_to_script; (2) qsub qsub_calculate_seed_profiles.py

# do relative import based on current working directory
# otherwise I have to install the package for relative import to work
current_wd = os.getenv('SGE_O_WORKDIR')
subpackage_folder_path = os.path.abspath(os.path.join(current_wd, '..'))
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import IO
import matchmaker
import type_conversions


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seed_folder", help="", type=str)
    parser.add_argument("--rna_bin_file", help="", type=str)
    parser.add_argument("--out_folder", help="", type=str)
    parser.add_argument("--inp_filename_template", help="", type=str)
    parser.add_argument("--out_filename_template", help="", type=str)
    parser.add_argument("--print_qstat", help="", type=str)
    parser.add_argument("--path_to_qstat", help="", type=str)



    parser.set_defaults(

        seed_folder='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/test_seeds',
        out_folder='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/test_profiles',
        inp_filename_template='test_seeds_101',
        out_filename_template='test_motifs_101',
        # seed_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/seeds/seeds_4-7_4-9_4-6_14-20/motifs_per_file_30k',
        # out_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/profiles/profiles_4-7_4-9_4-6_14-20/profiles_per_file_30k',
        # inp_filename_template='seeds_4-7_4-9_4-6_14-20_30k',
        # out_filename_template='profiles_4-7_4-9_4-6_14-20_30k',
        rna_bin_file='/wynton/home/goodarzi/khorms/pyteiser_root/data/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',


        c_command_path='/wynton/home/goodarzi/khorms/iTEISER/c_programs_new/calculate_profile_seed_june_18_mk',
        path_to_qstat='/opt/sge/bin/lx-amd64/qstat',
        print_qstat='y'

    )

    args = parser.parse_args()

    return args


def get_env_variables():
    working_home_dir = os.getcwd()
    local_scratch = os.getenv('TMPDIR')
    task_id = os.environ["SGE_TASK_ID"]
    job_id = os.environ["JOB_ID"]

    env_variables_dict = {"working_dir": working_home_dir,
                          "local_scratch": local_scratch,
                          "task_id": task_id,
                          "job_id": job_id}
    return env_variables_dict


def print_qstat_proc(env_variables_dict, args):
    #subprocess.call(qstat_command, shell=False)  # qstat is an executable itself, you don't need shell to run it
    # also subprocess doesn't want to find qstat if I provide the arguments as a string
    # for more detail, see either the comment of jfs here https://stackoverflow.com/questions/18962785/oserror-errno-2-no-such-file-or-directory-while-using-python-subprocess-in-dj
    # or this post https://stackoverflow.com/questions/4795190/passing-variables-to-a-subprocess-call

    subprocess.call([args.path_to_qstat, '-j', env_variables_dict["job_id"]], shell=False)  # qstat is an executable itself, you don't need shell to run it


def get_current_in_out_filenames(args, env_variables_dict):
    inp_filename_short = "%s_%s.bin" % (args.inp_filename_template, env_variables_dict["task_id"])
    out_filename_short = "%s_%s.bin" % (args.out_filename_template, env_variables_dict["task_id"])
    seeds_filename_full = os.path.join(args.seed_folder, inp_filename_short)
    profiles_filename_full = os.path.join(args.out_folder, out_filename_short)
    rna_bin_filename = args.rna_bin_file

    return seeds_filename_full, profiles_filename_full, rna_bin_filename



def calculate_write_profiles(n_motifs_list, n_seqs_list,
                            out_filename, do_print=False,
                             do_return = False):
    with open(out_filename, 'wb') as wf:
        if do_return:
            profiles_list = [0] * len(n_motifs_list)

        for i, motif in enumerate(n_motifs_list):
            current_profile, time_spent = matchmaker.calculate_profile_one_motif(motif, n_seqs_list)
            current_profile.compress()
            wf.write(current_profile.bytestring)

            if do_return:
                profiles_list[i] = current_profile.values

            if do_print:
                print("Motif number %d binds %d sequences. It took %.2f seconds"
                      % (i, current_profile.sum(), time_spent))
    if do_return:
        profiles_array = np.array(profiles_list)
        return profiles_array


def read_input_files(seeds_filename_full, rna_bin_filename):
    seqs_dict, seqs_order = IO.read_rna_bin_file(rna_bin_filename)
    w_motifs_list = IO.read_motif_file(seeds_filename_full)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]
    n_motifs_list = type_conversions.w_to_n_motifs_list(w_motifs_list)
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)

    return n_motifs_list, n_seqs_list



def main():
    args = handler()
    env_variables_dict = get_env_variables()
    seeds_filename_full, profiles_filename_full, rna_bin_filename = get_current_in_out_filenames(args, env_variables_dict)
    n_motifs_list, n_seqs_list = read_input_files(seeds_filename_full, rna_bin_filename)
    calculate_write_profiles(n_motifs_list, n_seqs_list,
                             profiles_filename_full, do_print=True)


    if args.print_qstat == 'y':
        print_qstat_proc(env_variables_dict, args)


if __name__ == "__main__":
    main()
