import numpy as np
import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import MI


def MI_get_pvalue_and_zscore(active_profile, discr_exp_profile,
                       current_MI, n_permutations):
    shuffled_MI_values = np.zeros(n_permutations, dtype=np.float64)
    shuffled_MI_sum = np.float64(0)
    shuffled_MI_sum_sqares = np.float64(0)

    for i in range(n_permutations):
        shuffled_expr = np.random.permutation(discr_exp_profile)
        ith_MI = MI.mut_info(active_profile, shuffled_expr)

        shuffled_MI_values[i] = ith_MI

        shuffled_MI_sum += ith_MI
        shuffled_MI_sum_sqares += ith_MI*ith_MI

    shuffled_MI_values.sort()

    if current_MI < shuffled_MI_values[0]:
        # shortcut: if current MI is less than the minimal permuted MI, exit
        value_undiv = n_permutations
    else:
        # go from right to left while the shuffled score is higher than the real one
        j = n_permutations - 1
        while (j >= 0) and (current_MI <= shuffled_MI_values[j]):
            j -= 1
        value_undiv = n_permutations - j - 1

    pvalue = value_undiv / float(n_permutations)

    # print(shuffled_MI_values)
    # print(current_MI)
    return pvalue

