import os
import argparse

import IO
import binarize_sequences
import preprocess_custom_expression_profile
import calculate_seed_profiles
import calculate_MI_profiles
import choose_significant_seeds_v3

def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rna_fastafile", help="fasta file with RNA sequences", type=str)
    parser.add_argument("--exp_values_file", help="expression values in a csv format", type=str)
    parser.add_argument("--anno_name_column", help="column name in exp_values file that contains annotations", type=str)
    parser.add_argument("--measur_column", help="column name in exp_values file that contains expression measurements", type=str)
    parser.add_argument("--seeds_file", help="file with seeds in binary format", type=str)
    parser.add_argument("--temp_folder", help="folder to write temporary files to", type=str)
    parser.add_argument("--out", help="output file", type=str)

    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)
    parser.add_argument("--min_occurences", help="minimal number of seed occurence in the transcriptome"
                                                 " for a seed to be considered", type=int)
    parser.add_argument("--n_permutations", help="number of permutations for the rank test for a seed", type=int)
    parser.add_argument("--max_pvalue", help="p-value threshold", type=float)
    parser.add_argument("--min_zscore", help="z-score threshold", type=float)
    parser.add_argument("--step_1_jump", help="step size at the 1st round of empirical threshold search", type=int)
    parser.add_argument("--step_2_min_interval", help="resolution of empirical threshold search", type=int)
    parser.add_argument("--step_1_min_fraction", help="minimal fraction of passing seeds for the 1st round of empirical threshold search", type=float)
    parser.add_argument("--step_2_min_fraction", help="minimal fraction of passing seeds for the 2nd round of empirical threshold search", type=float)
    parser.add_argument("--step_3_min_fraction", help="minimal fraction of passing seeds for the 3rd round of empirical threshold search", type=float)

    parser.add_argument("--are_input_seeds_degenerate", help="do input seeds contain degenerate nucleotides", type=bool)
    parser.add_argument("--indices_mode", help="compression in the index mode", type=bool)
    parser.add_argument("--index_bit_width", help="number of bits per one index when compressing", type=int)


    parser.set_defaults(
        rna_fastafile='/Users/student/Documents/hani/programs/pyteiser/data/tutorial_example_files/test_seqs.fa',
        exp_values_file='/Users/student/Documents/hani/programs/pyteiser/data/tutorial_example_files/test_expression_values.csv',
        anno_name_column='ensembl',
        measur_column='measurement',
        seeds_file='/Users/student/Documents/hani/programs/pyteiser/data/tutorial_example_files/test_seeds_20.bin',
        temp_folder='/Users/student/Documents/hani/programs/pyteiser/data/temp',
        out='/Users/student/Documents/hani/programs/pyteiser/data/test_output.txt',
        nbins=2,
        min_occurences=5,
        are_input_seeds_degenerate = False,
        n_permutations=1000,  # takes 1 second per 100 permutations, Hani's default number of permutations is 1*10^6
        max_pvalue=0.05,  # Hani's default threshold is 1*10^-7
        min_zscore=-1,

        step_1_jump=10,  # Hani's default jump is 200
        step_2_min_interval=2,

        step_1_min_fraction=0.8,
        step_2_min_fraction=0.8,
        step_3_min_fraction=0.9,

        indices_mode = True,
        index_bit_width = 24,
    )

    args = parser.parse_args()
    return args


def generate_filenames(args):
    # list filenames
    temp_folder = args.temp_folder
    inputs_folder = os.path.join(temp_folder, "inputs")
    interm_folder = os.path.join(temp_folder, "interm")
    rna_bin_filename = os.path.join(inputs_folder, "rna.bin")
    exp_mask_filename = os.path.join(inputs_folder, "exp_mask.bin")
    profiles_filename = os.path.join(interm_folder, "profiles.bin")
    MI_filename = os.path.join(interm_folder, "MI_values.bin")
    passed_seeds_filename = os.path.join(interm_folder, "passed_seeds.bin")
    passed_profiles_filename = os.path.join(interm_folder, "passed_profiles.bin")

    # make directories
    IO.create_folder(inputs_folder)
    IO.create_folder(interm_folder)

    filenames_dict = {
        "rna_bin" : rna_bin_filename,
        "exp_mask" : exp_mask_filename,
        "profiles" : profiles_filename,
        "MI" : MI_filename,
        "passed_seeds" : passed_seeds_filename,
        "passed_profiles" : passed_profiles_filename
    }

    return filenames_dict

def main():
    args = handler()
    filenames_dict = generate_filenames(args)

    binarize_sequences_args = ["--rna_fastafile", args.rna_fastafile, "--rna_bin_file", filenames_dict['rna_bin']]
    preprocess_custom_expression_profile_args = ["--rna_bin_file", filenames_dict['rna_bin'], "--exp_values_file",
        args.exp_values_file, "--exp_mask_file", filenames_dict['exp_mask'], "--anno_name_column",
        args.anno_name_column, "--measur_column", args.measur_column
    ]

    binarize_sequences.main(binarize_sequences_args)
    preprocess_custom_expression_profile.main(preprocess_custom_expression_profile_args)
    calculate_seed_profiles.non_sge_dependent_main(
                           args.seeds_file, filenames_dict['profiles'], filenames_dict['rna_bin'],
                           args.are_input_seeds_degenerate, args.indices_mode, args.index_bit_width)
    calculate_MI_profiles.non_sge_dependent_main(
                            filenames_dict['profiles'], filenames_dict['MI'],
                            filenames_dict['exp_mask'],
                            args.indices_mode,
                            args.nbins,
                            args.min_occurences,
                            do_print = True)
    choose_significant_seeds_v3.non_sge_dependent_main(args.seeds_file,
                           filenames_dict['profiles'],
                           filenames_dict['MI'],
                           filenames_dict['passed_seeds'],
                           filenames_dict['passed_profiles'],
                           filenames_dict['exp_mask'],
                           args.max_pvalue, args.min_zscore, args.n_permutations,
                           args.step_1_min_fraction, args.step_1_jump, args.step_2_min_interval,
                           args.step_3_min_fraction,
                           args.indices_mode,
                           args.index_bit_width,
                           do_print=True)


if __name__ == '__main__':
    main()

