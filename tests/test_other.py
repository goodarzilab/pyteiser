import os
import sys
import numpy as np
import numba

current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)

import pyteiser.numba_replacement_functions as numba_rf






def test_np_unique_return_counts():
    one_arr = np.array([1,2,3,3,2,1,2,2,2,1])
    two_arr = np.array([1,1,1,2,2,2,3,3,3,1])
    three_arr = one_arr + two_arr
    test_array_3D = np.stack((three_arr, one_arr, two_arr)).transpose()
    test_array_1D = one_arr.copy()

    exp_value_1D, exp_counts_1D = np.unique(test_array_1D, return_counts=True, axis=0)
    exp_value_3D, exp_counts_3D = np.unique(test_array_3D, return_counts=True, axis=0)

    # print(exp_value_1D, exp_counts_1D)
    # print(exp_value_3D, exp_counts_3D)

    test_value_1D, test_counts_1D = numba_rf.np_unique_return_counts(test_array_1D)
    test_value_3D, test_counts_3D = numba_rf.np_unique_return_counts(test_array_3D)


    assert(np.array_equal(exp_value_1D, test_value_1D))
    assert(np.array_equal(exp_counts_1D, test_counts_1D))
    assert(np.array_equal(exp_value_3D, test_value_3D))
    assert(np.array_equal(exp_counts_3D, test_counts_3D))

    # print(exp_value_1D, test_value_1D)
    # print(exp_counts_1D, test_counts_1D)
    # print(exp_value_3D, test_value_3D)
    # print(exp_counts_3D, test_counts_3D)


def main():
    test_np_unique_return_counts()


if __name__ == "__main__":
    main()



