import numpy as np
import hashlib
import numba

import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import glob_var


class w_motif:
    # w stands for wrapper. This is an external class that is used to interact with the environment, read, write etc
    # for fast operations, we use n_motif, an internal class that is compatible with numba

    def __init__(self, stem_length, loop_length):
        self.stem_length = np.uint8(stem_length)
        self.loop_length = np.uint8(loop_length)
        self.length = stem_length + loop_length
        self.linear_length = stem_length * 2 + loop_length
        self.sequence = np.ones(shape=self.length, dtype=np.uint8)
        self.structure = np.repeat(np.array([glob_var._stem, glob_var._loop], dtype=np.uint8),
                                   np.array([stem_length, loop_length], dtype=np.uint8))


    def print_sequence(self, return_string = False):
        string_to_print = ''.join([glob_var._char_to_nt_mapping[x] for x in self.sequence])
        if not return_string:
            print(string_to_print)
        else:
            return string_to_print


    def print_structure(self, return_string = False):
        string_to_print = ''.join([glob_var._char_to_struct_mapping[x] for x in self.structure])
        if not return_string:
            print(string_to_print)
        else:
            return string_to_print


    def print(self):
        self.print_sequence()
        self.print_structure()


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
        return motif_copy



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


spec_profile = [
    ('values', numba.uint8)
]


@numba.jitclass(spec_profile)
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

@numba.jitclass(spec_motif)
class n_motif:
    # n stands for numba
    # this class is used for fast calculations with numba
    # for input/output operations, use w_motif class

    def __init__(self, stem_length, loop_length, sequence, structure):
        self.stem_length = stem_length
        self.loop_length = loop_length
        self.length = stem_length + loop_length
        self.linear_length = stem_length * 2 + loop_length
        self.sequence = sequence
        self.structure = structure


spec_sequence = [
    ('length', numba.uint32),
    ('nts', numba.uint8[:])
]

@numba.jitclass(spec_sequence)
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


# Some functions that break numba compilation if I put them inside the class
@numba.jit(cache=True, nopython=True, nogil=True)
def copy_n_motif(motif):
    motif_copy = n_motif(motif.stem_length, motif.loop_length,
                         motif.sequence, motif.structure)
    return motif_copy