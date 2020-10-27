import os
import argparse

import IO
import binarize_sequences
import preprocess_custom_expression_profile
import calculate_seed_profiles

def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rna_fastafile", help="fasta file with RNA sequences", type=str)
    parser.add_argument("--exp_values_file", help="expression values in a csv format", type=str)
    parser.add_argument("--anno_name_column", help="column name in exp_values file that contains annotations", type=str)
    parser.add_argument("--measur_column", help="column name in exp_values file that contains expression measurements", type=str)
    parser.add_argument("--seeds_file", help="file with seeds in binary format", type=str)
    parser.add_argument("--temp_folder", help="folder to write temporary files to", type=str)
    parser.add_argument("--out", help="output file", type=str)

    parser.set_defaults(
        rna_fastafile='/Users/student/Documents/hani/programs/pyteiser/data/tutorial_example_files/test_seqs.fa',
        exp_values_file='/Users/student/Documents/hani/programs/pyteiser/data/tutorial_example_files/test_expression_values.csv',
        anno_name_column='ensembl',
        measur_column='measurement',
        seeds_file='/Users/student/Documents/hani/programs/pyteiser/data/tutorial_example_files/test_seeds.bin',
        temp_folder='/Users/student/Documents/hani/programs/pyteiser/data/temp',
        out='/Users/student/Documents/hani/programs/pyteiser/data/test_output.txt'
    )

    args = parser.parse_args()
    return args


def generate_filenames(args):
    # list filenames
    temp_folder = args.temp_folder
    inputs_folder = os.path.join(temp_folder, "inputs")
    rna_bin_filename = os.path.join(inputs_folder, "rna.bin")
    exp_mask_filename = os.path.join(inputs_folder, "exp_mask.bin")

    # make directories
    IO.create_folder(inputs_folder)

    filenames_dict = {
        "rna_bin" : rna_bin_filename,
        "exp_mask" : exp_mask_filename
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
    calculate_seed_profiles


if __name__ == '__main__':
    main()