import sys
import numpy as np
import numba
import scipy.stats
import math



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


def entropy(labels, how, base=None):
    value, counts = np.unique(labels, return_counts=True)
    if how == 1:
        _ = entropy1(counts)
    elif how == 2:
        _ = entropy2(counts, labels.shape[0])
    elif how == 3:
        _ = entropy3(counts)
    elif how == 4:
        _ = entropy4(counts, labels.shape[0])

    return _



# double entropy_empirical(std::map< std::vector<int> ,int > frequencies, int nb_samples) {
#       double e = 0;
#       for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
#             e -= iter->second * log((double)iter->second);
#       return log((double)nb_samples) + e/nb_samples;
# }

#def entropy_empirical():


# @numba.jit(cache=True, nopython=True, nogil=True)
def entropy1(counts, base=None):
  return scipy.stats.entropy(counts, base=base)



def entropy2(counts, total_number, base=None):
  """ Computes entropy of label distribution. """

  probs = np.divide(counts, total_number)
  ent = 0.
  base = math.e if base is None else base

  for i in probs:
    ent -= i * math.log(i, base)

  return ent

@numba.jit(cache=True, nopython=True, nogil=True)
def entropy3(counts, base=None):
  norm_counts = np.divide(counts, counts.sum())
  base = math.e if base is None else base
  return -(norm_counts * np.log(norm_counts)/np.log(base)).sum()





import timeit

def time_entropies():
    vect_to_discr_10k = np.random.normal(size=10000)
    discr_vect = discretize(vect_to_discr_10k, bins=5, disc="equalfreq")

    a1 = entropy(discr_vect, how=1)
    a2 = entropy(discr_vect, how=2)
    a3 = entropy(discr_vect, how=3)
    print(a1, a2, a3)

    time_entr_1 = timeit.timeit(lambda: entropy(discr_vect, how = 1), number=1000)
    time_entr_2 = timeit.timeit(lambda: entropy(discr_vect, how = 2), number=1000)
    time_entr_3 = timeit.timeit(lambda: entropy(discr_vect, how = 3), number=1000)
    print("Entropy 1: ", time_entr_1)
    print("Entropy 2: ", time_entr_2)
    print("Entropy 3: ", time_entr_3)





def main():
    time_entropies()





if __name__ == "__main__":
    main()

















