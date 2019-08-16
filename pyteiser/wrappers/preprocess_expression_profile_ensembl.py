import numpy as np
import pandas as pd
import argparse

import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import IO


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rna_bin_file", help="", type=str)
    parser.add_argument("--exp_df_ensembl", help="", type=str)
    parser.add_argument("--exp_mask_file", help="", type=str)

    parser.set_defaults(
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        exp_df_ensembl='/Users/student/Documents/hani/programs/pyteiser/data/expression_data/TARBP2_decay_t_score_ensembl.txt',
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin'
    )

    args = parser.parse_args()
    return args


def read_ensembl_df(args, return_meas_dict = True):
    exp_df_ensembl = pd.read_csv(args.exp_df_ensembl, sep='\t')
    exp_df_ensembl.index = exp_df_ensembl['ensembl']

    if return_meas_dict:
        measurements_dict_full = exp_df_ensembl.to_dict()
        measurements_dict = measurements_dict_full['measurement']

        return exp_df_ensembl, measurements_dict
    else:
        return exp_df_ensembl


def construct_mask_arrays(args):
    seqs_dict, seqs_order = IO.read_rna_bin_file(args.rna_bin_file)
    exp_df_ensembl, measurements_dict = read_ensembl_df(args)

    transcripts_measured_list = exp_df_ensembl['ensembl'].tolist()
    transcripts_measured_set = set(transcripts_measured_list)

    list_indices_occuring = [1 if x in transcripts_measured_set else 0 for x in seqs_order]
    list_measurement_values = [measurements_dict[x] if x in transcripts_measured_set else 0 for x in
                                seqs_order]

    array_indices_occuring = np.array(list_indices_occuring, dtype=np.bool)
    array_measurement_values = np.array(list_measurement_values, dtype=np.float32)

    return array_indices_occuring, array_measurement_values


def compress_write_mask_arrays(index_array, values_array, args):
    assert(index_array.shape == values_array.shape)
    length_uint32 = np.array([index_array.shape], dtype=np.uint32)
    length_bitstring = length_uint32.tobytes()
    index_array_bytes = index_array.tobytes()
    values_array_bytes = values_array.tobytes()
    full_bytes_string = length_bitstring + index_array_bytes + values_array_bytes

    with open(args.exp_mask_file, 'wb') as wb:
        wb.write(full_bytes_string)


def main():
    args = handler()
    array_indices_occuring, array_measurement_values = construct_mask_arrays(args)
    compress_write_mask_arrays(array_indices_occuring, array_measurement_values, args)



if __name__ == '__main__':
    main()
