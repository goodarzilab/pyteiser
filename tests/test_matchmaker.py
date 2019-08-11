# see how it's implemented here: https://www.kennethreitz.org/essays/repository-structure-and-python

import os
import sys
import numpy as np
sys.path.insert(0, os.path.abspath('..'))

import pyteiser.glob_var as glob_var
import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.matchmaker as matchmaker
import pyteiser.type_conversions as type_conversions


def test_matchmaker():
    # test matchmaking algorithms
    # 3 strings listed here contain instances of 3 matches that are also listed here

    test_motif_1 = structures.w_motif(4,6)
    test_motif_2 = structures.w_motif(4,6)
    test_motif_3 = structures.w_motif(4,6)
    test_motif_1.from_string("GNCANCNNUU")
    test_motif_2.from_string("AAUNNGNGNU")
    test_motif_3.from_string("NNACGNNCUU")
    test_motifs_list_w = [test_motif_1, test_motif_2, test_motif_3]
    test_motifs_list = type_conversions.w_to_n_motifs_list(test_motifs_list_w)

    test_string_1 = 'UUUUUUUGACAACAAUUTGTCUUUUU' # instance motif_1 at 7
    test_string_2 = "GGCAUCAGUUUUUUAAUGUGUGAUCAUUGGGUUCCCCCUUUUU" # instance motif_2 at 14
    test_string_3 = "AAUUAAAACCCCCCCAAACGCCCUUGUUUCCCACCACGGGCUUGUGGAAAAUUUUUU" # instances motif_3 at 15 and 33

    test_sequence_1 = structures.w_sequence(len(test_string_1))
    test_sequence_2 = structures.w_sequence(len(test_string_2))
    test_sequence_3 = structures.w_sequence(len(test_string_3))
    test_sequence_1.from_sequence(test_string_1)
    test_sequence_2.from_sequence(test_string_2)
    test_sequence_3.from_sequence(test_string_3)
    test_sequences_list_w = [test_sequence_1, test_sequence_2, test_sequence_3]
    test_sequences_list = type_conversions.w_to_n_sequences_list(test_sequences_list_w)

    boolean_matchmaker_desired = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=bool)
    boolean_matchmaker_res = np.zeros(shape=(3, 3), dtype=bool)
    indices_matchmaker_desired = [[7],[],[],[],[14],[],[],[],[15,33]]
    indices_matchmaker_res = []

    for i, mt in enumerate(test_motifs_list):
        for k, sq in enumerate(test_sequences_list):
            is_match = matchmaker.is_there_motif_instance(mt, sq)
            matching_indices = matchmaker.find_all_motif_instances(mt, sq)
            boolean_matchmaker_res[i,k] = is_match
            indices_matchmaker_res.append(matching_indices)


    assert(np.array_equal(boolean_matchmaker_res, boolean_matchmaker_desired))
    assert(indices_matchmaker_res == indices_matchmaker_desired)


if __name__ == "__main__":
    test_matchmaker()
