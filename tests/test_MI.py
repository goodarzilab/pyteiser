import os
import sys
import numpy as np
import numba
import argparse

current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)

import pyteiser.MI as MI




def test_discret_eq_freq():

    vect_to_discr_9 = np.array([0.5, 5.1, 5.2, 4.8, 9.9, 0.1, 9.7, 0.2, 10.3], dtype=np.float32)
    discr_expected_9_result = np.array([0,1,1,1,2,0,2,0,2], dtype=np.uint16) # 3 bins
    vect_to_discr_9_result = MI.discretize(vect_to_discr_9, bins=3, disc = "equalfreq")
    assert(np.array_equal(vect_to_discr_9_result, discr_expected_9_result))

    vect_to_discr_30 = np.array([-0.490, 1.761, -1.400, 0.411, -0.244, 0.177, 0.091,
                          -0.349, -0.554, 1.339, -0.094, 0.757, -0.469, -0.973, -1.192,
                          -0.831, -0.618, 0.335, 0.020, 0.406, 0.301, -1.721, -0.678,
                          -0.917, 1.498, -1.084, -0.152, -0.915, 0.094, -0.499], dtype=np.float32)
    discr_expected_30_result = np.array([2, 6, 0, 5, 3, 4, 4, 3, 2, 6, 3, 6, 2, 0, 0, 1, 1, 5,
                                         4, 5, 5, 0, 1, 1, 6, 0, 3, 1, 4, 2], dtype=np.uint16)  # 7 bins
    vect_to_discr_30_result = MI.discretize(vect_to_discr_30, bins=7, disc = "equalfreq")
    assert (np.array_equal(vect_to_discr_30_result, discr_expected_30_result))


def test_discret_eq_width():

    vect_to_discr_9 = np.array([0.5, 5.1, 5.2, 4.8, 9.9, 0.1, 9.7, 0.2, 10.3], dtype=np.float32)
    discr_expected_9_result = np.array([0,1,1,1,2,0,2,0,2], dtype=np.uint16) # 3 bins
    vect_to_discr_9_result = MI.discretize(vect_to_discr_9, bins=3, disc = "equalwidth")
    assert(np.array_equal(vect_to_discr_9_result, discr_expected_9_result))

    vect_to_discr_30 = np.array([-0.490, 1.761, -1.400, 0.411, -0.244, 0.177, 0.091,
                          -0.349, -0.554, 1.339, -0.094, 0.757, -0.469, -0.973, -1.192,
                          -0.831, -0.618, 0.335, 0.020, 0.406, 0.301, -1.721, -0.678,
                          -0.917, 1.498, -1.084, -0.152, -0.915, 0.094, -0.499], dtype=np.float32)
    discr_expected_30_result = np.array([2, 6, 0, 4, 2, 3, 3, 2, 2, 6, 3, 4, 2, 1, 1, 1, 2, 4, 3,
                                         4, 4, 0, 2, 1, 6, 1, 3, 1, 3, 2], dtype=np.uint16)  # 7 bins
    vect_to_discr_30_result = MI.discretize(vect_to_discr_30, bins=7, disc = "equalwidth")
    assert (np.array_equal(vect_to_discr_30_result, discr_expected_30_result))


def test_mutinf():
    one_arr = np.array([1, 2, 3, 3, 2, 1, 2, 2, 2, 1])
    two_arr = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3, 1])

    mi_test = MI.mut_info(one_arr, two_arr, with_numba=False)
    mi_expected = 0.28418101912817351
    #print(mi_test, mi_expected)
    assert(np.isclose(mi_test, mi_expected, atol=1e-16))

    cmi_test = MI.cond_mut_info(one_arr, two_arr, one_arr + two_arr)
    cmi_expected = 0.50219293007150134
    #print(cmi_test, cmi_expected)
    assert(np.isclose(cmi_test, cmi_expected, atol=1e-16))

    # an example from here: https://nlp.stanford.edu/IR-book/html/htmledition/mutual-information-1.html#mifeatsel2
    ut = np.repeat([0, 1, 0, 1], [774106, 27625, 141, 49])
    cc = np.repeat([0, 0, 1, 1], [774106, 27625, 141, 49])
    mi_test_base_2 = MI.mut_info(ut, cc, base=2, with_numba=False)
    mi_expected_base_2 = 0.0001105
    #print(mi_test_base_2, mi_expected_base_2)
    assert(np.isclose(mi_test_base_2, mi_expected_base_2, atol=1e-6))


def main():
    test_discret_eq_freq()
    test_discret_eq_width()
    test_mutinf()


if __name__ == "__main__":
    main()
