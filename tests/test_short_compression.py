import argparse
import numpy as np
import timeit
import hashlib
import struct
import numba

import os
import sys

# to make sure relative import works in order to import test data
current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)
os.chdir(package_home_path)


import pyteiser.glob_var as glob_var
import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.matchmaker as matchmaker



def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seedfile", type=str)
    parser.add_argument("--rna_fastafile", type=str)
    parser.add_argument("--outfile", type=str)
    parser.add_argument("--profiles_full_file", type=str)
    parser.add_argument("--compressed_profiles_file", type=str)
    parser.add_argument("--indices_mode", help="compression in the index mode", type=bool)

    parser.set_defaults(
        rna_fastafile='tests/data/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.txt',
        rna_bin_file='tests/data/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        profiles_full_file='tests/data/test_1_2_profiles_unique.bin',
        compressed_profiles_file='tests/data/compressed_by_indices_profiles.bin',
        indices_mode=False
    )

    args, unknown = parser.parse_known_args()

    return args


def run_test_bins_fasta(args):
    with open(args.rna_bin_file, 'rb') as rb:
        bitstring = rb.read()
        seq_objects_dict, seq_objects_order = IO.decompress_named_sequences(bitstring)
        full_string = IO.write_named_seq_to_fasta(seq_objects_dict, seq_objects_order)
    with open(args.rna_fastafile, 'r') as rf:
        full_fasta_string = rf.read()


    full_fasta_string_Us = full_fasta_string.replace('T','U').replace('ENSU','ENST')
    assert(len(full_string) == len(full_fasta_string_Us))

    assert(full_string == full_fasta_string_Us)



def run_test_compressing_decompressing_indices(args, decompressed_profiles_array):
    with open(args.compressed_profiles_file, 'wb') as wf:
        transcriptome_length = decompressed_profiles_array.shape[1]
        for i in range(decompressed_profiles_array.shape[0]):
            curr_profile = structures.w_profile(transcriptome_length)
            curr_profile.values = decompressed_profiles_array[i]
            curr_profile.compress_indices()
            wf.write(curr_profile.bytestring_indices)

    with open(args.compressed_profiles_file, 'rb') as rf:
        bitstring = rf.read()

    read_out_profiles_array = IO.decompress_profiles_indices(bitstring)
    assert (read_out_profiles_array == decompressed_profiles_array).all(), "decompression has changed the data!"


def profiles_wrapper(args):
    profiles_filename = args.profiles_full_file
    decompressed_profiles_array = IO.unpack_profiles_file(profiles_filename,
                                                          args.indices_mode,
                                                          do_print = True)
    run_test_compressing_decompressing_indices(args, decompressed_profiles_array)


def test_main():
    args = handler()
    run_test_bins_fasta(args)
    profiles_wrapper(args)


if __name__ == "__main__":
    test_main()

