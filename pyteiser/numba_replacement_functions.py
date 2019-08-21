import numpy as np
import numba
from numba import numpy_support


def np_unique_return_counts(inp_ar):
    ar = inp_ar.copy()
    orig_shape, orig_dtype = ar.shape, ar.dtype
    # reshape -1 says: numpy, figure it out yourself
    # see https://stackoverflow.com/questions/18691084/what-does-1-mean-in-numpy-reshape
    ar = ar.reshape(orig_shape[0], -1)

    aux = ar.copy()
    aux.sort(axis=0)

    # find where are there new values
    mask = np.empty(aux.shape[0], dtype=np.bool)
    mask[:1] = True
    mask[1:] = (aux[1:, :] != aux[:-1, :]).any(axis=1)

    indices_temp = np.nonzero(mask)[0]  # np.nonzero returns a tuple
    indices = np.zeros(indices_temp.shape[0] + 1, dtype=np.int)
    indices[:-1] = indices_temp
    indices[-1] = mask.size

    unique_ones = aux[mask]
    # I don't really know how does this reshaping work, I took
    # this line from the reshape_uniq function
    # from here: https://github.com/numpy/numpy/blob/v1.17.0/numpy/lib/arraysetops.py
    unique_ones = unique_ones.reshape(-1, *orig_shape[1:])

    counts = np.diff(indices)

    return unique_ones, counts


np_unique_return_counts(U)
np_unique_return_counts(one_arr)


def main():
    pass


main()