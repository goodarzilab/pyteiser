import argparse
import numpy as np
import timeit
import hashlib
import struct
import numba

import os
import sys

# to make sure relative import works
# for a detailed explanation, see test_matchmaker.py

current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)

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

    parser.set_defaults(
        rna_fastafile = '/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/narrow_down_transcripts_list/Gencode_v28_GTEx_expressed_transcripts_fasta/utr_3_fasta/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.txt',
        # rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin',
        rna_bin_file='/Users/student/Documents/hani/temp/temp_bins/test_bin.bin',
        #profiles_full_file='/Users/student/Documents/hani/programs/pyteiser/data/test_1_batch_tarbp2/profiles_4-7_4-9_4-6_14-20_30k_1.bin',
        profiles_full_file='/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique.bin',
        compressed_profiles_file='/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/compressed_by_indices_profiles.bin',
    )

    args = parser.parse_args()

    return args


def test_bins_fasta(args):
    with open(args.rna_bin_file, 'rb') as rb:
        bitstring = rb.read()
        seq_objects_dict, seq_objects_order = IO.decompress_named_sequences(bitstring)
        full_string = IO.write_named_seq_to_fasta(seq_objects_dict, seq_objects_order)
    with open(args.rna_fastafile, 'r') as rf:
        full_fasta_string = rf.read()

    with open("/Users/student/Documents/hani/temp/temp_fasta/1.txt",'w') as wf:
        wf.write(full_string)

    full_fasta_string_Us = full_fasta_string.replace('T','U').replace('ENSU','ENST')
    assert(len(full_string) == len(full_fasta_string_Us))
    #
    # print(full_string[0:200])
    # print(full_fasta_string_Us[0:200])
    assert(full_string == full_fasta_string_Us)

# def compress_2(pr, width = 20):
#     # compression that takes less space for sparse profiles
#     # instead of saving each element of an array as 1 byte, we save indices of all the elements that are True
#     # saving one index requires more than 16 bytes (since 16 bytes allow for 65536 values, and a transcriptome can be larger)
#     # however, using 32 bytes per index seems like a waste of space
#     # therefore, I'll be using 20 bytes per index; this allows for 1 mln values
#     #
#     # the format is:
#     # first, 4 bytes keep the length of the profile in uint32 format
#     # then, 4 bytes keep total number of indexes of True values (N)
#     # then, one array of length 20 * N keeps indices positions
#     # then, last 16 bytes keep MD5 checksum
#
#     true_indices = np.where(pr) # get indices
#     true_indices = true_indices[0] # for some reason np.where returns a tuple
#     true_indices = true_indices.astype(np.uint32) # convert to unsigned integers
#
#     # split each uint32 value into 4 uint8 values so that we can turn them into binary in a vectorized manner
#     N_indices = true_indices.shape[0]
#     binary_array = np.zeros((N_indices, 32), dtype=np.bool)
#     # iterate through all the indices and turn each into bin and then into 4 uint8 values
#     for i, value in enumerate(true_indices):
#         bit_string = struct.pack('I', value)
#         curr_uint8 = np.frombuffer(bit_string, dtype=np.uint8)
#         binary_array[i, :] = np.unpackbits(curr_uint8)
#
#     # remove extra unused bytes: shorten each value from 32 bits to the specified width (20 by default)
#     total_count_per_byte = binary_array.sum(axis=0)
#     assert (total_count_per_byte[width:] == 0).all() , "some indices are larger than the chosen width!"
#     shortened_binary_array = binary_array[:, 0:width]
#     flattened_binary_array = shortened_binary_array.flatten()
#
#     # make the total array of K*8 length so that we can compress it to bytes
#     total_number_of_values = flattened_binary_array.shape[0]
#     if total_number_of_values % 8 != 0:
#         new_number_of_values = ((total_number_of_values // 8) + 1) * 8
#         binary_bytes_array = np.zeros(new_number_of_values, dtype=np.bool)
#         binary_bytes_array[0:total_number_of_values] = flattened_binary_array
#     else:
#         binary_bytes_array = flattened_binary_array
#
#     length_uint32 = np.array([pr.shape[0]], dtype=np.uint32)
#     length_bitstring = length_uint32.tobytes()
#
#     N_indices_uint32 = np.array([N_indices], dtype=np.uint32)
#     N_indices_bitstring = N_indices_uint32.tobytes()
#
#     indices_packbits = np.packbits(binary_bytes_array)
#     indices_bitstring = indices_packbits.tobytes()
#
#     info_bitstring = length_bitstring + N_indices_bitstring + indices_bitstring
#
#     md5 = hashlib.md5()
#     md5.update(info_bitstring)
#     md5_checksum = md5.digest()
#     assert (md5.digest_size == 16)  # md5 checksum is always 16 bytes long, see wiki: https://en.wikipedia.org/wiki/MD5
#
#     full_bytestring = info_bitstring + md5_checksum
#     return full_bytestring
    # pr.bytestring = sequence_bytestring
    # pr.md5 = md5_checksum


def decompress_short_compression(bitstring):
    profiles_list = []
    total_length = len(bitstring)
    current_spot = 0
    counter = 0

    while current_spot < total_length:
        # get the length of the profile
        length_bitstring = bitstring[current_spot : current_spot + 4]
        profile_length_np = np.frombuffer(length_bitstring, dtype=np.uint32)
        length = profile_length_np[0]

        # get the number of indices (of True) of the profile
        N_indices_bitstring = bitstring[current_spot + 4 : current_spot + 8]
        N_indices_np = np.frombuffer(N_indices_bitstring, dtype=np.uint32)
        N_indices = N_indices_np[0]

        # get the number of bits used per index (compression width)
        width_bitstring = bitstring[current_spot + 8 : current_spot + 12]
        width_np = np.frombuffer(width_bitstring, dtype=np.uint32)
        width = width_np[0]

        # figure out how many bytes do we need to read out
        length_packed = N_indices * width

        if length_packed % 8 != 0:
            length_packed = (length_packed // 8) + 1
        else:
            length_packed = length_packed // 8

        # read out bitstring of the proper size
        values_bitstring = bitstring[current_spot + 12 : current_spot + 12 + length_packed]
        md5_bitstring = bitstring[current_spot + 12 + length_packed :
                                    current_spot + 12 + length_packed + 16]
        current_spot += 12 + length_packed + 16

        # convert bitsting to 32-bit arrays representing indices
        indices_packed_uint8 = np.frombuffer(values_bitstring, dtype=np.uint8)
        binary_bytes_array = np.unpackbits(indices_packed_uint8)
        binary_bytes_array = binary_bytes_array[0 : N_indices * width]
        reshaped_binary_array = binary_bytes_array.reshape(N_indices, width)
        full_binary_array = np.zeros((N_indices, 32), dtype=np.bool)
        full_binary_array[:, 0:width] = reshaped_binary_array

        # convert 32-bit arrays into a uint32 indices
        reshaped_full_binary_array = full_binary_array.flatten()
        reshaped_full_binary_string = np.packbits(reshaped_full_binary_array)
        true_indices = np.frombuffer(reshaped_full_binary_string, dtype=np.uint32)

        # create a new profile
        curr_profile = structures.w_profile(length)
        curr_profile.values[true_indices] = True
        curr_profile.compress_indices()
        assert (md5_bitstring == curr_profile.md5_indices)
        profiles_list.append(curr_profile.values)

        counter += 1
        # if counter % how_often_print == 0:
        #     if do_print:
        #         print("Decompressed profile number ", counter)

    profiles_array = np.array(profiles_list, dtype=np.bool)

    return profiles_array




    #
    #     values = np.unpackbits(values_packed_bits)
    #     values = values[0 : profile_length]
    #
    #     current_profile = structures.w_profile(profile_length)
    #     current_profile.values = values
    #     current_profile.compress()
    #
    #     assert (md5_bitstring == current_profile.md5)
    #
    #     profiles_list.append(current_profile.values)
    #
    #     counter += 1
    #     if counter % how_often_print == 0:
    #         if do_print:
    #             print("Decompressed profile number ", counter)
    #
    # profiles_array = np.array(profiles_list, dtype=np.bool)
    #
    # return profiles_array



def test_compressing_decompressing_indices(args, decompressed_profiles_array):
    with open(args.compressed_profiles_file, 'wb') as wf:
        transcriptome_length = decompressed_profiles_array.shape[1]
        for i in range(decompressed_profiles_array.shape[0]):
            curr_profile = structures.w_profile(transcriptome_length)
            curr_profile.values = decompressed_profiles_array[i]
            curr_profile.compress_indices()
            wf.write(curr_profile.bytestring_indices)

    with open(args.compressed_profiles_file, 'rb') as rf:
        bitstring = rf.read()

    read_out_profiles_array = decompress_short_compression(bitstring)
    assert (read_out_profiles_array == decompressed_profiles_array).all(), "decompression has changed the data!"


def profiles_wrapper(args):
    profiles_filename = args.profiles_full_file
    decompressed_profiles_array = IO.unpack_profiles_file(profiles_filename, do_print = True)
    current_array = decompressed_profiles_array[6, :]
    test_compressing_decompressing_indices(args, decompressed_profiles_array)




def main():
    args = handler()
    #test_bins_fasta(args)
    profiles_wrapper(args)


if __name__ == "__main__":
    main()

