import numpy as np
import argparse

import os
import sys


MASK_OUT_SEED_VALUE = np.float64(-1)


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task_mapping_file", help="", type=str)

    parser.add_argument("--profiles_folder", help="", type=str)
    parser.add_argument("--MI_values_folder", help="", type=str)
    parser.add_argument("--inp_filename_template", help="", type=str)
    parser.add_argument("--out_filename_template", help="", type=str)

    parser.add_argument("--print_qstat", help="", type=str)
    parser.add_argument("--path_to_qstat", help="", type=str)

    parser.add_argument("--rna_bin_file", help="", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)

    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)
    parser.add_argument("--min_occurences", help="minimal number of seed occurence in the transcriptome"
                                                 " for a seed to be considered", type=int)

    parser.add_argument("--calculate_with_numba", help="which MI calculation implementation to use: "
                                                       "numba-based or numba-free", type=str)


    parser.set_defaults(
        profiles_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/profiles/profiles_4-7_4-9_4-6_14-20/profiles_per_file_30k',
        MI_values_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/MI_values/MI_values_4-7_4-9_4-6_14-20/MI_values_per_file_30k',
        inp_filename_template='profiles_4-7_4-9_4-6_14-20_30k',
        out_filename_template='MI_values_4-7_4-9_4-6_14-20_30k',

        rna_bin_file='/wynton/home/goodarzi/khorms/pyteiser_root/data/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        exp_mask_file='/wynton/home/goodarzi/khorms/pyteiser_root/data/mask_files/TARBP2_decay_t_score_mask.bin',


        path_to_qstat='/opt/sge/bin/lx-amd64/qstat',
        print_qstat='y',

        nbins = 15,
        min_occurences = 5,

        calculate_with_numba='y',
    )

    args = parser.parse_args()

    return args


def import_modules():
    current_wd = os.getenv('SGE_O_WORKDIR')
    subpackage_folder_path = os.path.abspath(os.path.join(current_wd, '..'))
    if subpackage_folder_path not in sys.path:
        sys.path.append(subpackage_folder_path)

    global MI
    global IO
    global sge

    import MI
    import IO
    import sge


def get_current_in_out_filenames(args, env_variables_dict, mapping_dict):
    file_index_to_use =  mapping_dict[env_variables_dict["task_id"]]
    inp_filename_short = "%s_%s.bin" % (args.inp_filename_template, file_index_to_use)
    out_filename_short = "%s_%s.bin" % (args.out_filename_template, file_index_to_use)
    profiles_filename_full = os.path.join(args.profiles_folder, inp_filename_short)
    MI_values_filename_full = os.path.join(args.MI_values_folder, out_filename_short)
    rna_bin_filename = args.rna_bin_file

    return profiles_filename_full, MI_values_filename_full, rna_bin_filename


def calculate_MI_for_seeds(decompressed_profiles_array, index_array, discr_exp_profile,
                       min_occurences, calculate_with_numba, do_print=False):
    if calculate_with_numba == 'yes' or calculate_with_numba == 'y':
        with_numba = True
    else:
        with_numba = False
    MI_values_array = np.zeros(decompressed_profiles_array.shape[0], dtype=np.float32)

    for i, profile in enumerate(decompressed_profiles_array):
        active_profile = profile[index_array]

        if active_profile.sum() <= min_occurences:
            MI_values_array[i] = MASK_OUT_SEED_VALUE
            # print("The seed number %d binds only %d transcripts" % (i, active_profile.sum()))
            continue

        MI_values_array[i] = MI.mut_info(active_profile, discr_exp_profile, with_numba)

        if do_print:
            if i % 1000 == 0 and i > 0:
                print("Profile number %d has been calculated" % i)

    MI_values_array = np.array(MI_values_array, dtype=np.float64) # make sure all elements are of the same size
    return MI_values_array


def main():
    # I only import things if I run this script itself
    # do relative import based on current working directory
    # otherwise I have to install the package for relative import to work
    import_modules()

    args = handler()

    # get mapping of task ids to input files
    mapping_dict = sge.parse_task_mapping_file(args.task_mapping_file)
    # get the task id
    env_variables_dict = sge.get_env_variables()
    # get the names of input and output files
    profiles_filename_full, MI_values_filename_full, rna_bin_filename = get_current_in_out_filenames(args, env_variables_dict, mapping_dict)

    decompressed_profiles_array, index_array, values_array = IO.unpack_profiles_and_mask(profiles_filename_full,
                                                                                         args.exp_mask_file, do_print=True)

    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, args.nbins)

    MI_values_array = calculate_MI_for_seeds(decompressed_profiles_array, index_array, discr_exp_profile,
                                         args.min_occurences, args.calculate_with_numba, do_print = True)
    IO.write_MI_values(MI_values_array, args.nbins, MI_values_filename_full)

    if args.print_qstat == 'y':
        sge.print_qstat_proc(env_variables_dict, args.path_to_qstat)


if __name__ == "__main__":
    main()
