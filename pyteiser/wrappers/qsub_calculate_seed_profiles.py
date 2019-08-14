#$ -S /wynton/home/goodarzi/khorms/miniconda3/envs/pyteiser_env/bin/python
#$ -t 1-2
#$ -l h_rt=00:29:00
#$ -l mem_free=15G
#$ -l scratch=15G
#$ -e /wynton/scratch/khorms/logs/test_stderr.txt
#$ -o /wynton/scratch/khorms/logs/test_stdout.txt
#$ -q short.q
#$ -r yes

import numpy as np
import argparse
import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)


import glob_var
import structures
import IO
import matchmaker
import type_conversions


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seed_folder", help="", type=str)
    parser.add_argument("--fasta_file", help="", type=str)
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


def prepare_lists_for_calculations(args):
    seqs_dict, seqs_order = IO.read_rna_bin_file(args.rna_bin_file)
    w_motifs_list = IO.read_motif_file(args.seedfile)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]
    n_motifs_list = type_conversions.w_to_n_motifs_list(w_motifs_list)
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)
    return n_motifs_list, n_seqs_list



def main():
    args = handler()
    n_motifs_list, n_seqs_list = prepare_lists_for_calculations(args)

    matchmaker.calculate_profiles_list_motifs(n_motifs_list, n_seqs_list, do_print=True)


if __name__ == "__main__":
    main()
