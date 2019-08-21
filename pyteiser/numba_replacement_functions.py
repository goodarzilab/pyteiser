import numpy as np
import numba



def unpack_tuple(x):
    """ Unpacks one-element tuples for use as return values """
    if len(x) == 1:
        return x[0]
    else:
        return x


def reshape_uniq(uniq, orig_dtype, orig_shape, axis):
    uniq = uniq.view(orig_dtype)
    uniq = uniq.reshape(-1, *orig_shape[1:])
    uniq = np.swapaxes(uniq, 0, axis)
    return uniq


def unique1d_return_counts(ar):
    ar = np.asanyarray(ar).flatten()

    ar.sort()
    aux = ar
    mask = np.empty(aux.shape, dtype=np.bool_)
    mask[:1] = True
    mask[1:] = aux[1:] != aux[:-1]

    ret = (aux[mask],)
    idx = np.concatenate(np.nonzero(mask) + ([mask.size],))
    ret += (np.diff(idx),)
    return ret


def np_unique_return_counts(inp_ar):
    ar = inp_ar.copy()
    axis = 0
    ar = np.swapaxes(ar, axis, 0)
    orig_shape, orig_dtype = ar.shape, ar.dtype
    ar = ar.reshape(orig_shape[0], -1)
    ar = np.ascontiguousarray(ar)
    dtype = [('f{i}'.format(i=i), ar.dtype) for i in range(ar.shape[1])]
    consolidated = ar.view(dtype)

    output = unique1d_return_counts(consolidated)
    output = (reshape_uniq(output[0], orig_dtype, orig_shape, axis),) \
             + output[1:]

    return unpack_tuple(output)


def main():
    one_arr = np.array([1,2,3,3,2,1,2,2,2,1])
    two_arr = np.array([1,1,1,2,2,2,3,3,3,1])
    three_arr = one_arr + two_arr
    U = np.stack((three_arr, one_arr, two_arr)).transpose()
    print(U)

    a,b = np_unique_return_counts(U)

    print(a)
    print(b)


main()