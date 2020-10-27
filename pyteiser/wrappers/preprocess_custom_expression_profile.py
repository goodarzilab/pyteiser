import numpy as np
import pandas as pd
import argparse

import os
import sys

from .. import IO


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rna_bin_file", help="binarized sequence file", type=str)
    parser.add_argument("--exp_values_file", help="expression values in a csv format", type=str)
    parser.add_argument("--exp_mask_file", help="output file: indicates which sequences are present in the "
                                                "expression file and the expression values for these sequences", type=str)

    parser.add_argument("--anno_name_column", help="column name in exp_values file that contains annotations", type=str)
    parser.add_argument("--measur_column", help="column name in exp_values file that contains expression measurements", type=str)

    parser.set_defaults(
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/SNRNPA1_SE.hg38.fl250.bin',
        exp_values_file='/Users/student/Documents/hani/programs/pyteiser/data/expression_data/hg38_miso_se.txt',
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/SNRNPA1_PSI_mask.bin',

        anno_name_column='eid',
        measur_column='diff',
    )

    args = parser.parse_args()
    return args


def read_exp_values_file(args, return_meas_dict = True):
    exp_df = pd.read_csv(args.exp_values_file, sep='\t',
                         dtype = {args.anno_name_column : str})
    exp_df.index = exp_df[args.anno_name_column]

    if return_meas_dict:
        measurements_dict_full = exp_df.to_dict()
        measurements_dict = measurements_dict_full[args.measur_column]

        return exp_df, measurements_dict
    else:
        return exp_df


def construct_mask_arrays(args):
    seqs_dict, seqs_order = IO.read_rna_bin_file(args.rna_bin_file)
    exp_df, measurements_dict = read_exp_values_file(args)

    transcripts_measured_list = exp_df[args.anno_name_column].tolist()
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
