import numpy as np
import numba
from numba import numpy_support


@numba.jit(cache=True, nopython=True, nogil=True)
def np_unique_return_counts(inp_ar):
    ar = inp_ar.copy()
    orig_shape, orig_dtype = ar.shape, ar.dtype
    # reshape -1 says: numpy, figure it out yourself
    # see https://stackoverflow.com/questions/18691084/what-does-1-mean-in-numpy-reshape
    ar = ar.reshape(orig_shape[0], -1)

    aux = ar.copy()

    # sort by the columns, first to last
    # from here: https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
    # see J.J's answer
    for i in range(ar.shape[1], 0, -1):
        aux = aux[aux[:, i - 1].argsort(kind='mergesort')]

    # find where are there new values
    mask = np.empty(aux.shape[0], dtype=np.bool_)  # for some reason doesn't work if I do np.bool
    mask[:1] = True

    # any(axis) is not implemented in numba
    bool_complete = aux[1:, :] != aux[:-1, :]

    # iteration over 2D array is not implemented in numba
    for index in range(bool_complete.shape[0]):
        mask[1 + index:] = bool_complete[index].any()

    # mask[1:] = (aux[1:,:] != aux[:-1,:]).any(axis=1)

    indices_temp = np.nonzero(mask)[0]  # np.nonzero returns a tuple
    indices = np.zeros(indices_temp.shape[0] + 1, dtype=np.int64)
    indices[:-1] = indices_temp
    indices[-1] = mask.size

    unique_ones = aux[mask]
    # I don't really know how does this reshaping work, I took
    # this line from the reshape_uniq function
    # from here: https://github.com/numpy/numpy/blob/v1.17.0/numpy/lib/arraysetops.py
    unique_ones = unique_ones.reshape(-1, *orig_shape[1:])

    counts = np.diff(indices)

    return unique_ones, counts