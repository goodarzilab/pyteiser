import numpy as np
import hashlib

import glob_var

class s_motif:

    def __init__(self, stem_length, loop_length):
        self.stem_length = stem_length
        self.loop_length = loop_length
        self.length = stem_length + loop_length
        self.sequence = np.ones(shape=self.length, dtype=np.uint8)
        self.structure = np.repeat([glob_var._stem, glob_var._loop], [stem_length, loop_length], axis=None)


    def print_sequence(self):
        string_to_print = ''.join([glob_var._char_to_nt_mapping[x] for x in self.sequence])
        print(string_to_print)

    def print_structure(self):
        string_to_print = ''.join([glob_var._char_to_struct_mapping[x] for x in self.structure])
        print(string_to_print)

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
        assert(md5.digest_size == 16)# md5 checksum is always 16 bytes long, see wiki: https://en.wikipedia.org/wiki/MD5

        motif_bytestring = motif_info + md5_checksum
        self.bytestring = motif_bytestring
        self.md5 = md5_checksum


