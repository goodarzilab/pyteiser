import sys
import os
import numpy as np
import numba
import math

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)


@numba.jit(cache=True, nopython=True, nogil=True)
def discretize(inp_array, bins, noise_std = 0.000000001, new_seed = False):
    if not new_seed:
        np.random.seed(57)
    length = len(inp_array)
    to_discr = inp_array + np.random.normal(0, noise_std, length)

    # got the idea from here: https://stackoverflow.com/questions/39418380/histogram-with-equal-number-of-points-in-each-bin
    bins_for_discr = np.interp(np.linspace(0, length, bins + 1),
                               np.arange(length),
                               np.sort(to_discr))
    bins_for_discr[-1] += 1 # otherwize numpy creates one extra bin with only 1 point
    digitized = np.digitize(to_discr, bins_for_discr)
    digitized = digitized - 1

    return digitized


# to check resolution of the respective data types
# do np.finfo(np.float32).resolution: 1e-06
# or np.finfo(np.float64).resolution: 1e-06
@numba.jit(cache=True, nopython=True, nogil=True)
def entropy_empirical(counts, total_number, base=None):
    probs = np.divide(counts, total_number)
    ent = 0.
    base = math.e if base is None else base

    for i in probs:
        if i == 0: # np.isclose is not supported by numba
            continue
        ent -= i * math.log(i) / math.log(base)

    return ent


def mut_info(X, Y, x_bins, y_bins, base=None):
    c_xy = np.histogram2d(X, Y, [x_bins, y_bins])[0]
    flatten_c_xy = c_xy.flatten()
    c_x = c_xy.sum(axis=1)
    c_y = c_xy.sum(axis=0)

    Hyx = entropy_empirical(flatten_c_xy, flatten_c_xy.sum(), base=base)
    Hx = entropy_empirical(c_x, c_x.sum(), base=base)
    Hy = entropy_empirical(c_y, c_y.sum(), base=base)
    res = Hx + Hy - Hyx
    if res < 0:
        res = 0

    return res


@numba.jit(cache=True, nopython=True, nogil=True)
def discretize_exp_profile(index_array, values_array, nbins):
    active_values_array = values_array[index_array]
    quant_values_array = discretize(active_values_array, bins=nbins)
    return quant_values_array

# # input: 3 one-dimensional arrays
# @numba.jit(cache=True, nopython=True, nogil=True)
# def cond_mut_info(X, Y, Z, base=None):
#     U = np.stack((X, Z, Y)).transpose()
#     Hyzx = entropy(U, base)
#     Hzx = entropy(U[: , 0:2], base)
#     Hyz = entropy(U[: , 0:3], base)
#     Hz = entropy(Z, base)
#     Ires = Hyz - Hz - Hyzx + Hzx
#
#     return Ires




















