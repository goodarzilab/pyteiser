import sys
import os
import numpy as np
import numba
import math

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)


@numba.jit(cache=True, nopython=True, nogil=True)
def discretize(inp_array, bins, noise_std = 0.000000001, new_seed = False):
    if not new_seed:
        np.random.seed(57)
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
        if i == 0: # np.isclose is not supported by numba
            continue
        ent -= i * math.log(i) / math.log(base)

    return ent


@numba.jit(cache=True, nopython=True, nogil=True)
def histogram_1D(X, x_bins):
    x_bins = np.int64(x_bins) # for some reason it doesn't work with np.int32
    hist_array = np.zeros(x_bins)
    for x_val in X:
        x_coord = int(x_val)
        hist_array[x_coord] += 1
    return hist_array



@numba.jit(cache=True, nopython=True, nogil=True)
def histogram_2D(X, Y, x_bins, y_bins):
    x_bins = np.int64(x_bins) # for some reason it doesn't work with np.int32
    y_bins = np.int64(y_bins)
    hist_array = np.zeros((x_bins, y_bins))
    for x_val, y_val in zip(X,Y):
        x_coord = int(x_val)
        y_coord = int(y_val)
        hist_array[x_coord, y_coord] += 1
    return hist_array


@numba.jit(cache=True, nopython=True, nogil=True)
def histogram_3D(X, Y, Z, x_bins, y_bins, z_bins):
    x_bins = np.int64(x_bins) # for some reason it doesn't work with np.int32
    y_bins = np.int64(y_bins)
    z_bins = np.int64(z_bins)
    hist_array = np.zeros((x_bins, y_bins, z_bins))
    for x_val, y_val, z_val in zip(X,Y,Z):
        x_coord = int(x_val)
        y_coord = int(y_val)
        z_coord = int(z_val)
        hist_array[x_coord, y_coord, z_coord] += 1
    return hist_array


# in current implementation, all vallues in both vectors have to be integers between 0 and x/y_bins
# this solution increases the speed of histogram step (the slowest step) 10-fold, but it sacrifises some
# flexibility
@numba.jit(cache=True, nopython=True, nogil=True)
def mut_info(X, Y, x_bins, y_bins, base=None):
    #c_xy = np.histogram2d(X, Y, [x_bins, y_bins])[0]
    c_xy = histogram_2D(X, Y, x_bins, y_bins)
    flatten_c_xy = c_xy.flatten()
    c_x = c_xy.sum(axis=1)
    c_y = c_xy.sum(axis=0)

    Hyx = entropy_empirical(flatten_c_xy, flatten_c_xy.sum(), base=base)
    Hx = entropy_empirical(c_x, c_x.sum(), base=base)
    Hy = entropy_empirical(c_y, c_y.sum(), base=base)
    res = Hx + Hy - Hyx
    if res < 0:
        res = 0

    return res


# this is based on calculations in Infotheo R package source code
# link: https://cran.r-project.org/web/packages/infotheo/index.html
@numba.jit(cache=True, nopython=True, nogil=True)
def cond_mut_info(X, Y, Z, x_bins, y_bins, z_bins, base=None):
    c_yz = histogram_2D(Y, Z, y_bins, z_bins)
    flatten_c_yz = c_yz.flatten()
    H_yz = entropy_empirical(flatten_c_yz, flatten_c_yz.sum(), base=base)

    c_xz = histogram_2D(X, Z, x_bins, z_bins)
    flatten_c_xz = c_xz.flatten()
    H_xz = entropy_empirical(flatten_c_xz, flatten_c_xz.sum(), base=base)

    c_z = histogram_1D(Z, z_bins)
    flatten_c_z = c_z.flatten()
    H_z = entropy_empirical(flatten_c_z, flatten_c_z.sum(), base=base)

    c_xyz = histogram_3D(X, Y, Z, x_bins, y_bins, z_bins)
    flatten_c_xyz = c_xyz.flatten()
    H_xyz = entropy_empirical(flatten_c_xyz, flatten_c_xyz.sum(), base=base)

    Ires = H_yz - H_z - H_xyz + H_xz

    return Ires


@numba.jit(cache=True, nopython=True, nogil=True)
def discretize_exp_profile(index_array, values_array, nbins):
    active_values_array = values_array[index_array]
    quant_values_array = discretize(active_values_array, bins=nbins)
    return quant_values_array




















