import numpy as np
import hashlib
import numba
import struct
import sys

import glob_var


class w_motif:
    # w stands for wrapper. This is an external class that is used to interact with the environment, read, write etc
    # for fast operations, we use n_motif, an internal class that is compatible with numba

    def __init__(self, stem_length, loop_length):
        self.stem_length = np.uint8(stem_length)
        self.loop_length = np.uint8(loop_length)
        self.length = stem_length + loop_length
        self.sequence = np.ones(shape=self.length, dtype=np.uint8)
        self.structure = np.repeat(np.array([glob_var._stem, glob_var._loop], dtype=np.uint8),
                                   np.array([stem_length, loop_length], dtype=np.uint8))
        self.adjust_linear_length()


    def print_sequence(self, return_string = False):
        string_to_print = ''.join([glob_var._char_to_nt_mapping[x] for x in self.sequence])
        if not return_string:
            print(string_to_print)
        else:
            return string_to_print


    def print_linear_sequence(self, return_string = False):
        self.get_linear_sequence()
        string_to_print = ''.join([glob_var._char_to_nt_mapping[x] for x in self.linear_sequence])
        if not return_string:
            print(string_to_print)
        else:
            return string_to_print

    def print_linear_structure(self, return_string = False):
        linear_structure_array = [None] * self.linear_length

        left_index = 0
        right_index = left_index + self.linear_length - 1

        for i in range(self.length):
            if self.structure[i] == glob_var._stem:
                linear_structure_array[right_index] = '>'
                linear_structure_array[left_index] = '<'
                left_index += 1
                right_index -= 1
            else:
                linear_structure_array[left_index] = '.'
                left_index += 1

        string_to_print = ''.join(linear_structure_array)

        if not return_string:
            print(string_to_print)
        else:
            return string_to_print


    def encode_linear_structure(self):
        full_structure_array = np.ones(shape=self.linear_length, dtype=np.uint8)

        left_index = 0
        right_index = left_index + self.linear_length - 1

        for i in range(self.length):
            if self.structure[i] == glob_var._stem:
                full_structure_array[right_index] = glob_var._right_stem
                full_structure_array[left_index] = glob_var._left_stem
                left_index += 1
                right_index -= 1
            else:
                full_structure_array[left_index] = glob_var._loop
                left_index += 1

        self.full_structure_encoding = full_structure_array


    def print_structure(self, return_string = False):
        string_to_print = ''.join([glob_var._char_to_struct_mapping[x] for x in self.structure])
        if not return_string:
            print(string_to_print)
        else:
            return string_to_print


    def print(self):
        self.print_sequence()
        self.print_structure()


    def print_linear(self):
        self.print_linear_sequence()
        self.print_linear_structure()


    def eq(self, other):
        return ((self.stem_length == other.stem_length) and
                (self.loop_length == other.loop_length) and
                np.array_equal(self.sequence, other.sequence) and
                np.array_equal(self.structure, other.structure))


    def from_string(self, string):
        upper_string = string.upper()
        try:
            assert(len(upper_string) == self.length)
        except AssertionError:
            print("Error: the length of the string is %d and the length of the sequence object is %d" %
                  (len(upper_string), self.length))
            sys.exit(1) # see here https://stackoverflow.com/questions/438894/how-do-i-stop-a-program-when-an-exception-is-raised-in-python

        for ind, let in enumerate(upper_string):
            try:
                current_nt = glob_var._nt_to_char_mapping[let]
            except KeyError:
                print("Error: this string contains a nucleotide that is not ACGT/U. More specifically, it's ", let)
                print("Exiting!")
                sys.exit(1)
            np.put(self.sequence, ind, current_nt)


    def compress(self):
        # byte string representation of the motif
        # first, 2 bytes keep stem_length and loop_length
        # then, two consecutive arrays hold sequence and structure
        # then, last 16 bytes keep MD5 checksum

        characteristic_numbers = np.array([self.stem_length, self.loop_length], dtype=np.uint8)
        characteristic_numbers_bitstring = characteristic_numbers.tobytes()

        sequence_bytes = self.sequence.tobytes()
        structure_bytes = self.structure.tobytes()

        motif_info = characteristic_numbers_bitstring + sequence_bytes + structure_bytes

        md5 = hashlib.md5()
        md5.update(motif_info)
        md5_checksum = md5.digest()
        assert (md5.digest_size == 16)  # md5 checksum is always 16 bytes long, see wiki: https://en.wikipedia.org/wiki/MD5

        motif_bytestring = motif_info + md5_checksum

        self.bytestring = motif_bytestring
        self.md5 = md5_checksum


    def copy(self):
        motif_copy = w_motif(self.stem_length, self.loop_length)
        motif_copy.sequence = self.sequence
        motif_copy.structure = self.structure
        motif_copy.adjust_linear_length()
        return motif_copy


    def adjust_linear_length(self):
        stem_count = np.sum(self.structure == glob_var._stem)
        loop_count = np.sum(self.structure == glob_var._loop)
        self.linear_length = 2 * stem_count + loop_count


    def change_structure_position(self, position, st_type):
        if st_type == glob_var._loop:
            self.structure[position] = glob_var._loop
        elif st_type == glob_var._stem:
            self.structure[position] = glob_var._stem
        else:
            print("Inappropriate secondary structure value!")
            sys.exit(1)
        self.adjust_linear_length()


    def get_linear_sequence(self):
        self.linear_sequence = np.zeros(self.linear_length, dtype=np.uint8)

        left_index = 0
        right_index = left_index + self.linear_length - 1

        for i in range(self.length):
            current_nt = self.sequence[left_index]
            complementary_nt = glob_var._complementary_deg_nt_dict[current_nt]

            self.linear_sequence[left_index] = self.sequence[left_index]

            if self.structure[i] == glob_var._stem:
                self.linear_sequence[right_index] = complementary_nt
                left_index += 1
                right_index -= 1
            else:
                left_index += 1







class w_sequence:
    # w stands for wrapper. This is an external class that is used to interact with the environment, read, write etc
    # for fast operations, we use n_motif, an internal class that is compatible with numba

    def __init__(self, length):
        self.length = np.uint32(length)
        self.nts = np.zeros(shape=self.length, dtype=np.uint8)

    def from_sequence(self, string):
        upper_string = string.upper()
        try:
            assert(len(upper_string) == self.length)
        except AssertionError:
            print("Error: the length of the string is %d and the length of the sequence object is %d" %
                  (len(upper_string), self.length))
            sys.exit(1) # see here https://stackoverflow.com/questions/438894/how-do-i-stop-a-program-when-an-exception-is-raised-in-python

        for ind, let in enumerate(upper_string):
            try:
                current_nt = glob_var._nt_to_char_mapping[let]
            except KeyError:
                print("Error: this string contains a nucleotide that is not ACGT/U. More specifically, it's ", let)
                print("Exiting!")
                sys.exit(1)
            np.put(self.nts, ind, current_nt)

    def print(self, return_string = False):
        letters_array = np.vectorize(glob_var._char_to_nt_mapping.get)(self.nts)
        string_to_print = "".join(letters_array)
        if return_string:
            return string_to_print
        else:
            print(string_to_print)


    def print_sequence(self, beginning = 0, end = 0, return_string = False):
        if end == 0:
            end = self.length
        string_to_print = ''.join([glob_var._char_to_nt_mapping[x] for x in self.nts[beginning : end]])
        if not return_string:
            print(string_to_print)
        else:
            return string_to_print


    def compress(self):

        # byte string representation of the sequence
        # first, 4 bytes keep the length of the uint32 format
        # then, one array (uint8 format) holds the nts array
        # then, last 16 bytes keep MD5 checksum

        length_uint32 = np.array([self.length], dtype=np.uint32)
        length_bitstring = length_uint32.tobytes()

        nts_bytes = self.nts.tobytes()

        sequence_info = length_bitstring + nts_bytes

        md5 = hashlib.md5()
        md5.update(sequence_info)
        md5_checksum = md5.digest()
        assert (md5.digest_size == 16)  # md5 checksum is always 16 bytes long, see wiki: https://en.wikipedia.org/wiki/MD5

        sequence_bytestring = sequence_info + md5_checksum

        self.bytestring = sequence_bytestring
        self.md5 = md5_checksum


class w_profile:
    # w stands for wrapper. This is an external class that is used to interact with the environment, read, write etc
    # for fast operations, we use n_motif, an internal class that is compatible with numba

    def __init__(self, n_sequences):
        self.values = np.zeros(n_sequences, dtype=np.bool)

    def sum(self):
        return self.values.sum()


    def compress(self):
        # byte string representation of the sequence
        # first, 4 bytes keep the length of the profile in uint32 format
        # then, one array (uint8 format) holds the compressed (to bits) values array
        # then, last 16 bytes keep MD5 checksum

        values_packbits = np.packbits(self.values)

        # we keep the length of the actual array cause when we pack it
        # we fill the rest of the last byte with zeros, so we can recover the length
        # of the packed array from length of unpacked array but not vice versa
        length_uint32 = np.array([self.values.shape[0]], dtype=np.uint32)
        length_bitstring = length_uint32.tobytes()
        values_bytes = values_packbits.tobytes()
        sequence_info = length_bitstring + values_bytes

        md5 = hashlib.md5()
        md5.update(sequence_info)
        md5_checksum = md5.digest()
        assert (md5.digest_size == 16)  # md5 checksum is always 16 bytes long, see wiki: https://en.wikipedia.org/wiki/MD5

        sequence_bytestring = sequence_info + md5_checksum

        self.bytestring = sequence_bytestring
        self.md5 = md5_checksum


    def compress_indices(self, width = 24):
        # compression that takes less space for sparse profiles
        # instead of saving each element of an array as 1 byte, we save indices of all the elements that are True
        # saving one index requires more than 16 bytes (since 16 bytes allow for 65536 values, and a transcriptome can be larger)
        # however, using 32 bytes per index seems like a waste of space
        # therefore, I'll be using 20 bytes per index; this allows for 1 mln values
        #
        # the format is:
        # first, 4 bytes keep the length of the profile in uint32 format
        # then, 4 bytes keep total number of indexes of True values (N)
        # then, one array of length 20 * N keeps indices positions
        # then, last 16 bytes keep MD5 checksum

        true_indices = np.where(self.values) # get indices
        true_indices = true_indices[0] # for some reason np.where returns a tuple
        true_indices = true_indices.astype(np.uint32) # convert to unsigned integers

        # split each uint32 value into 4 uint8 values so that we can turn them into binary in a vectorized manner
        N_indices = true_indices.shape[0]
        binary_array = np.zeros((N_indices, 32), dtype=np.bool)
        # iterate through all the indices and turn each into bin and then into 4 uint8 values

        for i, value in enumerate(true_indices):
            bit_string = struct.pack('I', value)
            curr_uint8 = np.frombuffer(bit_string, dtype=np.uint8)
            binary_array[i, :] = np.unpackbits(curr_uint8)

        # remove extra unused bytes: shorten each value from 32 bits to the specified width (20 by default)
        total_count_per_byte = binary_array.sum(axis=0)
        assert (total_count_per_byte[width:] == 0).all() , "some indices are larger than the chosen width!"
        shortened_binary_array = binary_array[:, 0:width]
        flattened_binary_array = shortened_binary_array.flatten()

        # make the total array of K*8 length so that we can compress it to bytes
        total_number_of_values = flattened_binary_array.shape[0]
        if total_number_of_values % 8 != 0:
            new_number_of_values = ((total_number_of_values // 8) + 1) * 8
            binary_bytes_array = np.zeros(new_number_of_values, dtype=np.bool)
            binary_bytes_array[0:total_number_of_values] = flattened_binary_array
        else:
            binary_bytes_array = flattened_binary_array

        length_uint32 = np.array([self.values.shape[0]], dtype=np.uint32)
        length_bitstring = length_uint32.tobytes()

        N_indices_uint32 = np.array([N_indices], dtype=np.uint32)
        N_indices_bitstring = N_indices_uint32.tobytes()

        width_uint32 = np.array([width], dtype=np.uint32)
        width_bitstring = width_uint32.tobytes()

        indices_packbits = np.packbits(binary_bytes_array)
        indices_bitstring = indices_packbits.tobytes()

        info_bitstring = length_bitstring + N_indices_bitstring + width_bitstring + indices_bitstring

        md5 = hashlib.md5()
        md5.update(info_bitstring)
        md5_checksum = md5.digest()
        assert (md5.digest_size == 16)  # md5 checksum is always 16 bytes long, see wiki: https://en.wikipedia.org/wiki/MD5
        full_bytestring = info_bitstring + md5_checksum

        self.bytestring_indices = full_bytestring
        self.md5_indices = md5_checksum


spec_profile = [
    ('values', numba.uint8)
]


@numba.experimental.jitclass(spec_profile)
class n_profile:
    # n stands for numba
    # this class is used for fast calculations with numba
    # for input/output operations, use w_motif class

    def __init__(self, values):
        self.values = values



spec_motif = [
    ('stem_length', numba.uint8),
    ('loop_length', numba.uint8),
    ('length', numba.uint8),
    ('linear_length', numba.uint8),
    ('sequence', numba.uint8[:]),
    ('structure', numba.uint8[:])
]

@numba.experimental.jitclass(spec_motif)
class n_motif:
    # n stands for numba
    # this class is used for fast calculations with numba
    # for input/output operations, use w_motif class

    def __init__(self, stem_length, loop_length, sequence, structure):
        self.stem_length = stem_length
        self.loop_length = loop_length
        self.length = stem_length + loop_length
        self.sequence = sequence
        self.structure = structure
        self.adjust_linear_length()

    def adjust_linear_length(self):
        stem_count = np.sum(self.structure == glob_var._stem)
        loop_count = np.sum(self.structure == glob_var._loop)
        self.linear_length = 2 * stem_count + loop_count

    def change_structure_position(self, position, st_type):
        if st_type == glob_var._loop:
            self.structure[position] = glob_var._loop
        elif st_type == glob_var._stem:
            self.structure[position] = glob_var._stem
        self.adjust_linear_length()


spec_sequence = [
    ('length', numba.uint32),
    ('nts', numba.uint8[:])
]

@numba.experimental.jitclass(spec_sequence)
class n_sequence:
    # n stands for numba
    # this class is used for fast calculations with numba
    # for input/output operations, use w_motif class

    def __init__(self, length, nts):
        self.length = length
        self.nts = nts


    def is_paired(self, left_index, right_index):
        base1 = self.nts[left_index]
        base2 = self.nts[right_index]

        if (base1 == glob_var._U and base2 == glob_var._A) or \
                (base1 == glob_var._C and base2 == glob_var._G) or \
                (base1 == glob_var._G and base2 == glob_var._C) or \
                (base1 == glob_var._A and base2 == glob_var._U) or \
                (base1 == glob_var._U and base2 == glob_var._G) or \
                (base1 == glob_var._G and base2 == glob_var._U):
            return True
        return False


    def nt_is_a(self, ind, nt):
        if nt == glob_var._N:
            return True
        else:
            base = self.nts[ind]
            if nt == base:
                return True
            return False


    def nt_is_a_degenerate(self, ind, nt_degen):
        if nt_degen == glob_var._N:
            return True
        elif nt_degen == glob_var._Y:
            if self.nts[ind] == glob_var._U or self.nts[ind] == glob_var._C:
                return True
            return False
        elif nt_degen == glob_var._R:
            if self.nts[ind] == glob_var._A or self.nts[ind] == glob_var._G:
                return True
            return False
        elif nt_degen == glob_var._K:
            if self.nts[ind] == glob_var._U or self.nts[ind] == glob_var._G:
                return True
            return False
        elif nt_degen == glob_var._M:
            if self.nts[ind] == glob_var._A or self.nts[ind] == glob_var._C:
                return True
            return False
        elif nt_degen == glob_var._S:
            if self.nts[ind] == glob_var._G or self.nts[ind] == glob_var._C:
                return True
            return False
        elif nt_degen == glob_var._W:
            if self.nts[ind] == glob_var._A or self.nts[ind] == glob_var._U:
                return True
            return False
        elif nt_degen == glob_var._B:
            if self.nts[ind] == glob_var._G or self.nts[ind] == glob_var._U or self.nts[ind] == glob_var._C:
                return True
            return False
        elif nt_degen == glob_var._D:
            if self.nts[ind] == glob_var._G or self.nts[ind] == glob_var._A or self.nts[ind] == glob_var._U:
                return True
            return False
        elif nt_degen == glob_var._H:
            if self.nts[ind] == glob_var._A or self.nts[ind] == glob_var._C or self.nts[ind] == glob_var._U:
                return True
            return False
        elif nt_degen == glob_var._V:
            if self.nts[ind] == glob_var._G or self.nts[ind] == glob_var._C or self.nts[ind] == glob_var._A:
                return True
            return False
        else:
            if nt_degen == self.nts[ind]:
                return True
            return False


# Get a deep copy of an n_motif
# I can't figure out how to make this function a part of n_motif class without numba crashing
# I can deep copy numpy arrays, like sequence and structure
# but I can't deepcopy features like stem_length and loop_length, because for small integers (less than 256)
# Python actually keeps only 1 copy of those in memory and they don't get copied
# see answer by Yu Hao here: https://stackoverflow.com/questions/31621997/why-deepcopy-of-list-of-integers-returns-the-same-integers-in-memory
#
# For copying arrays, I have two options: (1) I can use deepcopy, it's safer but it's not numba-compatible
# Alternatively, I can use np.copy(). The array copy created by np.copy doesn't change when the original
# array gets changed. See examples in the source code:
# https://github.com/numpy/numpy/blob/v1.17.0/numpy/lib/function_base.py#L745-L790
def copy_n_motif(motif):
    copy_sequence = np.copy(motif.sequence)  # alternatively, use copy.deepcopy
    copy_structure = np.copy(motif.structure)  # alternatively, use copy.deepcopy

    motif_copy = n_motif(motif.stem_length, motif.loop_length, copy_sequence, copy_structure)

    return motif_copy
