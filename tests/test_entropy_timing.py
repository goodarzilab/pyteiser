import scipy.stats
import math
import numpy as np
import numba
import os
import sys

# to make sure relative import works
# for a detailed explanation, see test_matchmaker.py

current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)

import pyteiser.MI as MI

# Here, I want to test what would be the fastest way to calculate entropy in Python
# There is a stackoverflow discussion of this problem here: https://stackoverflow.com/questions/15450192/fastest-way-to-compute-entropy-in-python
# A good set of tests is provided in Jarad's answer; however, I want to wrap the same functions to Numba and therefore
# Jarad's tests are not sufficient for me. Here, I test the same functions but wrapped in numba


def entropy(labels, how, base=None):
    value, counts = np.unique(labels, return_counts=True)
    if how == "scipy":
        res = entropy_scipy(counts)
    elif how == "math":
        res = entropy_math(counts, labels.shape[0])
    elif how == "numpy":
        res = entropy_numpy(counts)
    else:
        print("Choose a proper method of entropy calculation!")
        sys.exit(1)

    return res


# numba doesn't want to work with scipy
def entropy_scipy(counts, base=None):
  return scipy.stats.entropy(counts, base=base)


@numba.jit(cache=True, nopython=True, nogil=True)
def entropy_math(counts, total_number, base=None):
  probs = np.divide(counts, total_number)
  ent = np.float64(0)
  base = math.e if base is None else base

  for i in probs:
    ent -= i * math.log(i) / math.log(base)

  return ent


@numba.jit(cache=True, nopython=True, nogil=True)
def entropy_numpy(counts, base=None):
  norm_counts = np.divide(counts, counts.sum())
  base = math.e if base is None else base
  return -(norm_counts * np.log(norm_counts)/np.log(base)).sum()


import timeit

def time_entropies():
    vect_to_discr_10k = np.random.normal(size=10000)
    discr_vect = MI.discretize(vect_to_discr_10k, bins=10, disc="equalfreq")

    e1 = entropy(discr_vect, how="scipy")
    e2 = entropy(discr_vect, how="math")
    e3 = entropy(discr_vect, how="numpy")

    assert(np.isclose(e1, e2, atol=1e-16))
    assert (np.isclose(e1, e3, atol=1e-16))
    # print(e1, e2, e3)

    time_entr_1 = timeit.timeit(lambda: entropy(discr_vect, how = "scipy"), number=1000)
    time_entr_2 = timeit.timeit(lambda: entropy(discr_vect, how = "math"), number=1000)
    time_entr_3 = timeit.timeit(lambda: entropy(discr_vect, how = "numpy"), number=1000)
    print("Entropy calculation with scipy takes : ", time_entr_1)
    print("Entropy calculation with math takes : ", time_entr_2)
    print("Entropy calculation with numpy takes : ", time_entr_3)


# example timing
# Entropy calculation with scipy takes :  0.2700598379999999
# Entropy calculation with math takes :  0.20465292200000018
# Entropy calculation with numpy takes :  0.22385192200000015

# multiple tests show that numpy-based and math-based calculations take approximately the same amount of time,
# whereas scipy is slightly slower. Also math-based is usually a tiny bit faster than numpy but not by much


def main():
    time_entropies()

if __name__ == "__main__":
    main()