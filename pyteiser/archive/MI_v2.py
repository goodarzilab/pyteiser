import sys
import os
import numpy as np
import numba
import math

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import numba_replacement_functions as numba_rf


# The MI calculation implementation based on np.unique runs very slow
# when applied to long arrays (40k or longer). I have replaced this implementation
# with another one that is based on np.histogram2d and I have archived the 
# np.unique implementation


@numba.jit(cache=True, nopython=True, nogil=True)
def discretize(inp_array, bins, noise_std = 0.000000001):
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
    ent -= i * math.log(i) / math.log(base)

  return ent


# numba doesn't support np.unique with return_counts argument
# therefore, np.unique-based implementation of entropy calculation is hard to speed up
# I have re-implemented np.unique in a numba-compatible way in the numba_replacement_functions file
# entropy_no_numba is a valid implementation which is not supported by numba
# however, in certain scenarios it can be much faster so I leave it in here
def entropy_no_numba(labels, base=None):
    value, counts = np.unique(labels, return_counts=True, axis=0)
    res = entropy_empirical(counts, len(labels), base)
    return res


@numba.jit(cache=True, nopython=True, nogil=True)
def entropy(labels, base=None):
    value, counts = numba_rf.np_unique_return_counts(labels)
    res = entropy_empirical(counts, len(labels), base)
    return res


# input: 2 one-dimensional arrays
@numba.jit(cache=True, nopython=True, nogil=True)
def mut_info_with_numba(X, Y, base=None):
    U = np.stack((X, Y)).transpose()
    Hyx = entropy(U, base)
    Hx = entropy(X, base)
    Hy = entropy(Y, base)
    res = Hx + Hy - Hyx
    if res < 0:
        res = 0

    return res


# input: 2 one-dimensional arrays
def mut_info_no_numba(X, Y, base=None):
    U = np.stack((X, Y)).transpose()
    Hyx = entropy_no_numba(U, base)
    Hx = entropy_no_numba(X, base)
    Hy = entropy_no_numba(Y, base)
    res = Hx + Hy - Hyx
    if res < 0:
        res = 0

    return res


# for some arrays, MI calculations without using my ad-hoc Numba implementation actually works faster
# therefore, it will be optional if the user wants to use numba here or not
def mut_info(X, Y, with_numba=True, base=None):
    if with_numba:
        return mut_info_with_numba(X, Y, base)
    else:
        return mut_info_no_numba(X, Y, base)


# input: 3 one-dimensional arrays
@numba.jit(cache=True, nopython=True, nogil=True)
def cond_mut_info(X, Y, Z, base=None):
    U = np.stack((X, Z, Y)).transpose()
    Hyzx = entropy(U, base)
    Hzx = entropy(U[: , 0:2], base)
    Hyz = entropy(U[: , 0:3], base)
    Hz = entropy(Z, base)
    Ires = Hyz - Hz - Hyzx + Hzx

    return Ires


def discretize_exp_profile(index_array, values_array, nbins):
    active_values_array = values_array[index_array]
    quant_values_array = discretize(active_values_array, bins=nbins)
    return quant_values_array



















