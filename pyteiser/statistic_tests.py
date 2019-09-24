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
                    n_permutations, max_pvalue, n_samples,
                   fraction_retain, min_fraction_passed,
                   do_print = False):
    total_number_passed = 0

    for j in range(n_samples):
        full_indices_array = np.arange(active_profile.shape[0])
        how_many_keep = int(fraction_retain * active_profile.shape[0])
        subsampl_index_array = np.random.choice(full_indices_array, size=how_many_keep, replace=False)
        curr_profile = active_profile[subsampl_index_array]
        curr_exp_profile = discr_exp_profile[subsampl_index_array]
        curr_MI = MI.mut_info(curr_profile, curr_exp_profile, x_bins=2, y_bins=nbins)
        pvalue, z_score = MI_get_pvalue_and_zscore(curr_profile, discr_exp_profile, nbins,
                                                   curr_MI, n_permutations)
        if do_print:
            print("Iteration %d. p-value: %.5f; max_pvalue: %.5f, z-score: %.2f" % (j, pvalue, max_pvalue, z_score))
        if pvalue < max_pvalue:
            total_number_passed += 1

    fraction_passed = total_number_passed / float(n_samples)
    if do_print:
        print("%.2f subsamples passed the test; required fraction is %.2f" % (fraction_passed, min_fraction_passed))
    if fraction_passed >= min_fraction_passed:
        if do_print:
            print("Passed robustness test")
        return True
    else:
        if do_print:
            print("Did not pass robustness test")
        return False
