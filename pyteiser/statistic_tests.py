import numpy as np
import os
import sys
import math
import numba

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)


import MI


#@numba.jit(cache=True, nopython=True, nogil=True)
def MI_get_pvalue_and_zscore(active_profile, discr_exp_profile,
                        nbins, current_MI, n_permutations):
    shuffled_MI_values = np.zeros(n_permutations, dtype=np.float64)

    for i in range(n_permutations):
        shuffled_expr = np.random.permutation(discr_exp_profile)
        ith_MI = MI.mut_info(active_profile, shuffled_expr, x_bins=2, y_bins=nbins)

        shuffled_MI_values[i] = ith_MI

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
    z_score = (current_MI - np.mean(shuffled_MI_values)) / np.std(shuffled_MI_values)

    # print(shuffled_MI_values)
    # print(current_MI)
    return pvalue, z_score


def jackknife_test(active_profile, discr_exp_profile, nbins,
                    n_permutations,
                    n_samples, fraction_retain, min_fraction_passed):
    total_number_passed = 0

    for i in n_samples:
        subsampl_index_array = np.subsample
        curr_profile = active_profile[subsampl_index_array]
        curr_exp_profile = discr_exp_profile[subsampl_index_array]
        curr_MI = MI.mut_info(curr_profile, curr_exp_profile, x_bins=2, y_bins=nbins)
    # jn = 10
    # jn_t = 6
    # jn_f = 3
    # pass = 0
    # do jn iterations:
        # subsample 1-1/jn_f = 2/3 of all the data
        # for subsample: calculate MI
        # for subsample: calculate p-value
        # if pvalue < max_pvalue:
            # pass += 1
    # if pass > jn_t: passed
    pass