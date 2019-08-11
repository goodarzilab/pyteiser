import numpy as np
import hashlib
import sys

from numba import jitclass
from numba import uint8


import glob_var


spec_motif = [
    ('stem_length', uint8),
    ('loop_length', uint8),
    ('length', uint8),
    ('linear_length', uint8),
    ('sequence', uint8[:]),
    ('structure', uint8[:])
]

@jitclass(spec_motif)
class s_motif:

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


class s_sequence:

    def __init__(self, length):
        self.length = length
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


    def print_sequence(self, beginning = 0, end = 0, return_string = False):
        if end == 0:
            end = self.length
        string_to_print = ''.join([glob_var._char_to_nt_mapping[x] for x in self.nts[beginning : end]])
        if not return_string:
            print(string_to_print)
        else:
            return string_to_print
