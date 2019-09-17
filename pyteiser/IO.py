import numpy as np
import os
import sys
import hashlib

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import glob_var
import structures



def decompress_motifs_from_bitstring(bitstring):
    motifs_list = []
    total_length = len(bitstring)
    current_spot = 0

    while current_spot < total_length:
        stem_length = bitstring[current_spot]
        loop_length = bitstring[current_spot + 1]
        full_length = stem_length + loop_length


        curr_sequence = np.frombuffer(bitstring[current_spot + 2 : current_spot + 2 + full_length], dtype=np.uint8)
        curr_structure = np.frombuffer(bitstring[current_spot + 2 + full_length :
                                        current_spot + 2 + 2*full_length], dtype=np.uint8)
        md5_checksum = bitstring[current_spot + 2 + 2*full_length: current_spot + 2 + 2*full_length + 16]

        current_motif = structures.w_motif(stem_length, loop_length)
        current_motif.sequence = curr_sequence
        current_motif.structure = curr_structure
        current_motif.compress()

        # current_motif.print_sequence()
        # current_motif.print_structure()

        assert(md5_checksum == current_motif.md5)

        motifs_list.append(current_motif)
        current_spot += 2 + 2*full_length + 16

    return motifs_list

def read_rna_bin_file(inp_file):
    with open(inp_file, 'rb') as rb:
        bitstring = rb.read()
        seq_objects_dict, seq_objects_order = decompress_named_sequences(bitstring)
    return seq_objects_dict, seq_objects_order


def read_motif_file(inp_file):
    with open(inp_file, 'rb') as rf:
        full_bitstring = rf.read()
        motifs_list = decompress_motifs_from_bitstring(full_bitstring)
    return motifs_list


def read_fasta(infile, do_print = False, how_often_print = 1000):
    tr_dict_loc = {}
    seqs_order = []
    with open(infile, 'r') as f:
        split_string = f.read().split('>')
        for ind, entry in enumerate(split_string):
            if entry == '':
                continue
            seq_start = entry.find('\n')
            annotation = entry[:seq_start]
            sequence_string = entry[seq_start + 1:].replace('\n', '')
            current_sequence = structures.w_sequence(len(sequence_string))
            current_sequence.from_sequence(sequence_string)
            tr_dict_loc[annotation] = current_sequence
            seqs_order.append(annotation)
            if ind % how_often_print == 0:
                if do_print:
                    print("Read sequence number ", ind)

    return tr_dict_loc, seqs_order


def compress_named_sequences(seq_objects_dict, seqs_order,
                             do_print=False, how_often_print=1000):

    seq_batch_byte_list = []

    for ind, name in enumerate(seqs_order):
        current_byte_string = b''
        full_length_name = name.ljust(glob_var.MAX_SEQ_NAME_LENGTH)
        name_in_bytes = full_length_name.encode('utf-8')
        assert(len(name_in_bytes) == glob_var.MAX_SEQ_NAME_LENGTH)
        current_byte_string += name_in_bytes

        seq_objects_dict[name].compress()

        current_byte_string += seq_objects_dict[name].bytestring
        seq_batch_byte_list.append(current_byte_string)
        if ind % how_often_print == 0:
            if do_print:
                print("Compressed sequence number ", ind)
    seq_batch_byte_string = b''.join(seq_batch_byte_list)

    return seq_batch_byte_string


def decompress_named_sequences(bitstring,
                               do_print=False, how_often_print=10000):
    seq_objects_dict = {}
    seq_objects_order = []

    total_length = len(bitstring)
    current_spot = 0
    counter = 0

    while current_spot < total_length:
        name_in_bytes = bitstring[current_spot : current_spot + glob_var.MAX_SEQ_NAME_LENGTH]
        length_bitstring = bitstring[current_spot + glob_var.MAX_SEQ_NAME_LENGTH :
                                        current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4]
        seq_length_np = np.frombuffer(length_bitstring, dtype=np.uint32)
        seq_length = seq_length_np[0]

        sequence_bitstring = bitstring[current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4 :
                                        current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4 + seq_length]
        md5_bitstring = bitstring[current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4 + seq_length :
                                        current_spot + glob_var.MAX_SEQ_NAME_LENGTH + 4 + seq_length + 16]

        current_spot += glob_var.MAX_SEQ_NAME_LENGTH + 4 + seq_length + 16

        full_length_name = name_in_bytes.decode('utf-8')
        name = full_length_name.rstrip()

        current_sequence = structures.w_sequence(seq_length)
        current_sequence.nts = np.frombuffer(sequence_bitstring, dtype = np.uint8)
        current_sequence.compress()


        assert (md5_bitstring == current_sequence.md5)

        seq_objects_dict[name] = current_sequence
        seq_objects_order.append(name)

        counter += 1
        if counter % how_often_print == 0:
            if do_print:
                print("Compressed sequence number ", counter)

    return seq_objects_dict, seq_objects_order


def write_named_seq_to_fasta(seq_objects_dict, seq_objects_order):
    strings_list = []
    for ind, name in enumerate(seq_objects_order):
        seq_string = seq_objects_dict[name].print_sequence(return_string = True)
        strings_list.append(">%s\n%s\n" % (name, seq_string))
    full_string = ''.join(strings_list)
    return full_string


def decompress_profiles(bitstring,
                        do_print=False, how_often_print=10000):
    profiles_list = []
    total_length = len(bitstring)
    current_spot = 0
    counter = 0

    while current_spot < total_length:
        # get the length of the profile
        length_bitstring = bitstring[current_spot : current_spot + 4]
        profile_length_np = np.frombuffer(length_bitstring, dtype=np.uint32)
        profile_length = profile_length_np[0]

        # figure out how long is the profile packed into bits
        # if profile length // 8 > 0, it will take one additional byte
        if profile_length % 8 != 0:
            length_packed = (profile_length // 8) + 1
        else:
            length_packed = profile_length // 8

        values_bitstring = bitstring[current_spot + 4 : current_spot + 4 + length_packed]
        md5_bitstring = bitstring[current_spot + 4 + length_packed :
                                    current_spot + 4 + length_packed + 16]

        current_spot += 4 + length_packed + 16

        values_packed_bits = np.frombuffer(values_bitstring, dtype=np.uint8)
        values = np.unpackbits(values_packed_bits)
        values = values[0 : profile_length]

        current_profile = structures.w_profile(profile_length)
        current_profile.values = values
        current_profile.compress()

        assert (md5_bitstring == current_profile.md5)

        profiles_list.append(current_profile.values)

        counter += 1
        if counter % how_often_print == 0:
            if do_print:
                print("Decompressed profile number ", counter)

    profiles_array = np.array(profiles_list, dtype=np.bool)

    return profiles_array


def decompress_exp_mask_file(bitstring):
    length_bitstring = bitstring[0 : 4]
    mask_length_np = np.frombuffer(length_bitstring, dtype=np.uint32)
    mask_length = mask_length_np[0]

    index_bitstring = bitstring[4 : 4 + mask_length] # dtype = np.bool
    values_bitstring = bitstring[4 + mask_length : 4 + mask_length + mask_length*4] # dtype=np.float32

    index_array = np.frombuffer(index_bitstring, dtype=np.bool)
    values_array = np.frombuffer(values_bitstring, dtype=np.float32)

    return index_array, values_array


def unpack_profiles_and_mask(profiles_bin_file, exp_mask_file, do_print=False):
    with open(profiles_bin_file, 'rb') as rf:
        bitstring = rf.read()
        decompressed_profiles_array = decompress_profiles(bitstring)
        if do_print:
            print("%d profiles have been loaded" % len(decompressed_profiles_array))

    with open(exp_mask_file, 'rb') as rf:
        bitstring = rf.read()
        index_array, values_array = decompress_exp_mask_file(bitstring)
        if do_print:
            print("Expression values are provided for %d out of %d transcripts in the reference transcriptome" %
                  (index_array.sum(), index_array.shape[0]))

    try:
        assert (decompressed_profiles_array[0].shape[0] == index_array.shape[0])
    except AssertionError:
        print("Error: occurence profiles were calculated for some other reference transcriptome. The length of the "
              "profiles is %d and the length of the transcriptome provided is %d" %
              (decompressed_profiles_array[0].shape[0], index_array.shape[0]))
        sys.exit(1)

    return decompressed_profiles_array, index_array, values_array


def write_MI_values(MI_values_array, nbins, MI_values_file):
    length_uint32 = np.array([MI_values_array.shape[0]], dtype=np.uint32)
    length_bitstring = length_uint32.tobytes()
    nbins_uint32 = np.array([nbins], dtype=np.uint32)
    nbins_bitstring = nbins_uint32.tobytes()
    values_bytes = MI_values_array.tobytes()
    MI_values_bytestring = length_bitstring + nbins_bitstring + values_bytes

    with open(MI_values_file, 'wb') as wf:
        wf.write(MI_values_bytestring)


def decompres_MI_values(bitstring):
    length_nbins_bitstring = bitstring[0: 8]
    MI_array_length_nbins_np = np.frombuffer(length_nbins_bitstring, dtype=np.uint32)
    MI_array_length = MI_array_length_nbins_np[0]
    nbins = MI_array_length_nbins_np[1]
    MI_array_bitstring = bitstring[8 : 8 + MI_array_length * 8] # np.float64 takes 8 bytes
    # to check it, you could run print(np.dtype(np.float32).itemsize)
    MI_array = np.frombuffer(MI_array_bitstring, dtype=np.float64)
    return MI_array, nbins


def read_MI_values(MI_values_file):
    with open(MI_values_file, 'rb') as rf:
        bitstring = rf.read()
    MI_values_array, nbins = decompres_MI_values(bitstring)
    return MI_values_array, nbins


# def write_seed_significancy_threshold(last_positive_seed, threshold_file):
#     threshold_bytes = np.uint32(last_positive_seed).tobytes()
#     md5 = hashlib.md5()
#     md5.update(threshold_bytes)
#     md5_checksum = md5.digest()
#     byte_string = threshold_bytes + md5_checksum
#     with open(threshold_file, 'wb') as wf:
#         wf.write(byte_string)
#
#
# def read_seed_significancy_threshold(threshold_file):
#     with open(threshold_file, 'rb') as rf:
#         bitstring = rf.read()
#     threshold_bytes = bitstring[0: 4]
#     md5_checksum_saved = bitstring[4:]
#     threshold_value = np.frombuffer(threshold_bytes, dtype=np.uint32)
#     threshold_bytes_recode = np.uint32(threshold_value).tobytes()
#     md5 = hashlib.md5()
#     md5.update(threshold_bytes_recode)
#     md5_checksum = md5.digest()
#     assert(md5_checksum == md5_checksum_saved)
#     return threshold_value[0]
#
#
# def decompres_seed_threshold(bitstring):
#     length_nbins_bitstring = bitstring[0: 8]
#     MI_array_length_nbins_np = np.frombuffer(length_nbins_bitstring, dtype=np.uint32)
#     MI_array_length = MI_array_length_nbins_np[0]
#     nbins = MI_array_length_nbins_np[1]
#     MI_array_bitstring = bitstring[8 : 8 + MI_array_length * 8] # np.float64 takes 8 bytes
#     # to check it, you could run print(np.dtype(np.float32).itemsize)
#     MI_array = np.frombuffer(MI_array_bitstring, dtype=np.float64)
#     return MI_array, nbins

def read_seed_pass_individual_file(inp_filename):
    with open(inp_filename, 'rb') as rf:
        bitstring = rf.read()
    length_bytes = bitstring[0: 4]
    length_value = np.frombuffer(length_bytes, dtype=np.uint32)
    if length_value == 0:
        return []

    motifs_bitstring = bitstring[4: ]
    motifs_list = decompress_motifs_from_bitstring(motifs_bitstring)
    assert(len(motifs_list) == length_value)
    return motifs_list


