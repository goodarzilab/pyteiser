import numpy as np



def discret_eq_freq(inp_array, nbins):
    N = inp_array.shape[0]

    spl = np.zeros(nbins, dtype=np.float32) # split points
    res = np.zeros(N, dtype=np.uint16) # discretized vector

    sorted_column = np.sort(inp_array)

    # find split points

    freq = N // nbins
    mod = N % nbins
    splitpoint = freq - 1

    for i in range(nbins - 1):
        if mod > 0:
            spl[i] = sorted_column[splitpoint + 1]
            mod -= 1
        else:
            spl[i] = sorted_column[splitpoint]
        splitpoint += freq

    EPSILON = np.float32(0.01)  # constant to add to the end split point
    spl[nbins - 1] = sorted_column[N - 1] + EPSILON

    # identify the bin corresponding to each element of an array

    for s in range(N):
        bin = -1
        k = 0
        while (bin == -1) and (k < nbins):
            if inp_array[s] <= spl[k]:
                bin = k
            res[s] = bin
            k += 1

    return res


def main():
    # no nans - remove them beforehand
    vect_to_discr = [0.5, 5.1, 5.2, 4.8, 9.9, 0.1, 9.7, 0.2, 10.3]
    discr_expected_result = [1, 2, 2, 2, 3, 1, 3, 1, 3]
    discr_expected_result = [x-1 for x in discr_expected_result]
    discr_expected_result_np = np.array(discr_expected_result, dtype=np.uint16)

    array_to_discr = np.array(vect_to_discr, dtype=np.float32)
    array_result = discret_eq_freq(array_to_discr, nbins=3)
    assert(np.array_equal(discr_expected_result_np, array_result))

main()

















