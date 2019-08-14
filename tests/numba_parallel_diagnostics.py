import os
import sys

# to make sure relative import works
# for a detailed explanation, see test_matchmaker.py

current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)


import numpy as np
import numba


import pyteiser.glob_var as glob_var
import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.matchmaker as matchmaker
import pyteiser.type_conversions as type_conversions
import pyteiser.wrappers.calculate_seed_profiles as calculate_seed_profiles



def create_single_pair(stem = 4, loop = 7,
                      motif_str = "NGCAUNGNANN",
                      seq_str = "UGCAUUGUAUGUGUG"):
    test_motif = structures.w_motif(stem, loop)
    test_motif.from_string(motif_str)
    n_test_motif = type_conversions.w_to_n_motif(test_motif)

    test_sequence = structures.w_sequence(len(seq_str))
    test_sequence.from_sequence(seq_str)
    n_test_sequence = type_conversions.w_to_n_sequence(test_sequence)

    return n_test_motif, n_test_sequence




def diagnoze_single_pair_function():
    n_test_motif, n_test_sequence = create_single_pair()

    # Dizagnosing is_there_motif_instance
    # _ = matchmaker.is_there_motif_instance(n_test_motif, n_test_sequence)
    # matchmaker.is_there_motif_instance.parallel_diagnostics(level=4)
    # is_there_motif_instance can't be optimized by numba parallel

    # Diagnosing match_motif_seq
    # _ = matchmaker.match_motif_seq(n_test_motif, n_test_sequence, 0)
    # matchmaker.match_motif_seq.parallel_diagnostics(level=4)
    # match_motif_seq can't be optimized by numba parallel





if __name__ == "__main__":
    diagnoze_single_pair_function()