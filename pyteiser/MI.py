import sys
import numpy as np
import numba




# This function is written based on discretize function from infotheo package for R: https://cran.r-project.org/web/packages/infotheo/index.html
# The original core fucntions - discEF and discEW - are in the discretize.cpp file is in the package's source code
# the three key differences are:
# (1) this implementation assumes the array doesn't contain any nans - you should remove them beforehand
# (2) this implementation only works with 1D arrays

#disc is the name of the discretization method to be used :"equalfreq" or "equalwidth" (default : "equalfreq")

def discretize(inp_array, bins, disc = "equalfreq"):
    inp_array_32 = np.array(inp_array, dtype=np.float32)
    unique_elements = np.unique(inp_array_32).shape[0]
    try:
        assert (unique_elements > 1)
    except AssertionError:
        print("All the values in the array are equal!")
        sys.exit(1)

    if np.isnan(inp_array_32).any():
        print("The input array contains NaNs!")
        sys.exit(1)

    if disc == "equalfreq":
        res = discret_eq_freq(inp_array, nbins = bins)
    elif disc == "equalwidth":
        res = discret_eq_width(inp_array, nbins = bins)
    else:
        print("Choose an appropriate discretization method!")
        sys.exit(1)

    return res


# the original discEF function has a bug that makes the bins size slighly unequal:
# the increment of splitpoint function happens outside of if-else loop. If the condition is True, the next item
# of "col" array is incremented by (freq) and splitpoint variable is incremented for (freq-1) only. This leads to
# unequal size of resulting bins.
# In the current implementation this bug is fixed

@numba.jit(cache=True, nopython=True, nogil=True)
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

            splitpoint += freq + 1
        else:
            spl[i] = sorted_column[splitpoint]

            splitpoint += freq
        #splitpoint += freq

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


@numba.jit(cache=True, nopython=True, nogil=True)
def discret_eq_width(inp_array, nbins):
    N = inp_array.shape[0]

    res = np.zeros(N, dtype=np.uint16)  # discretized vector

    max = np.max(inp_array)
    min = np.min(inp_array)
    binsize = (max - min) / nbins

    for s in range(N):
        b = 0
        while not ((min + b * binsize <= inp_array[s]) and (inp_array[s] < min + (b + 1) * binsize)):
            b += 1

        if b == nbins:
            b = nbins - 1

        res[s] = b

    return res




# double entropy_empirical(std::map< std::vector<int> ,int > frequencies, int nb_samples) {
#       double e = 0;
#       for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
#             e -= iter->second * log((double)iter->second);
#       return log((double)nb_samples) + e/nb_samples;
# }

#def entropy_empirical():







def main():
    pass





if __name__ == "__main__":
    main()

















