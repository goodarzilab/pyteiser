import numpy as np
import timeit
import numba

_bin_U = np.uint8(1)
_bin_C = np.uint8(2)
_bin_G = np.uint8(4)
_bin_A = np.uint8(8)
_bin_N = np.uint8(15)


_bin_Y = np.uint8(3) # UC
_bin_R = np.uint8(12) # AG
_bin_K = np.uint8(5) # UG
_bin_M = np.uint8(10) # AC
_bin_S = np.uint8(6) # GC
_bin_W = np.uint8(9) # AU

_bin_B = np.uint8(7) # GUC
_bin_D = np.uint8(13) # GAU
_bin_H = np.uint8(11) # ACU
_bin_V = np.uint8(14) # GCA


def create_arrays_5_letters_only():
    real_let_only = [_bin_U, _bin_C, _bin_G, _bin_A]
    real_let_N = [_bin_U, _bin_C, _bin_G, _bin_A, _bin_N]
    all_deg_let = [_bin_U, _bin_C, _bin_G, _bin_A, _bin_N,
                   _bin_Y, _bin_R, _bin_K, _bin_M, _bin_S,
                   _bin_W, _bin_B, _bin_D, _bin_H, _bin_V]

    real_letters_list = np.ones(shape=120, dtype=np.uint8)
    deg_letters_list_N = np.ones(shape=120, dtype=np.uint8)
    deg_letters_list_full = np.ones(shape=120, dtype=np.uint8)

    for i, q in enumerate(real_let_N):
        for k in range(24):
            deg_letters_list_N[i * 24 + k] = q

    for i, q in enumerate(all_deg_let):
        for k in range(8):
            deg_letters_list_full[i * 8 + k] = q

    counter = 0
    for k in range(30):
        for q in real_let_only:
            real_letters_list[counter] = q
            counter += 1

    return real_letters_list, deg_letters_list_N, deg_letters_list_full


@numba.jit(cache=True, nopython=True, nogil=True)
def is_L_simple(L_pat, L_seq):
    if L_pat == _bin_N:
        return True
    else:
        if L_pat == L_seq:
            return True
        return False


@numba.jit(cache=True, nopython=True, nogil=True)
def is_L_deg(L_pat, L_seq):
    if L_pat == _bin_N:
        return True
    elif L_pat == _bin_Y:
        if L_seq == _bin_U or L_seq == _bin_C:
            return True
        return False
    elif L_pat == _bin_R:
        if L_seq == _bin_A or L_seq == _bin_G:
            return True
        return False
    elif L_pat == _bin_K:
        if L_seq == _bin_U or L_seq == _bin_G:
            return True
        return False
    elif L_pat == _bin_M:
        if L_seq == _bin_A or L_seq == _bin_C:
            return True
        return False
    elif L_pat == _bin_S:
        if L_seq == _bin_G or L_seq == _bin_C:
            return True
        return False
    elif L_pat == _bin_W:
        if L_seq == _bin_A or L_seq == _bin_U:
            return True
        return False
    elif L_pat == _bin_B:
        if L_seq == _bin_G or L_seq == _bin_U or L_seq == _bin_C:
            return True
        return False
    elif L_pat == _bin_D:
        if L_seq == _bin_G or L_seq == _bin_A or L_seq == _bin_U:
            return True
        return False
    elif L_pat == _bin_H:
        if L_seq == _bin_A or L_seq == _bin_C or L_seq == _bin_U:
            return True
        return False
    elif L_pat == _bin_V:
        if L_seq == _bin_G or L_seq == _bin_C or L_seq == _bin_A:
            return True
        return False
    else:
        if L_pat == L_seq:
            return True
        return False


@numba.jit(cache=True, nopython=True, nogil=True)
def bitwize_comp(i,k):
    return np.bitwise_and(i, k)


def compare_L_simple_pattern(pat_list, seq_list, do_return=False):
    counter = 0
    for i, k in zip(pat_list, seq_list):
        if is_L_simple(i, k):
            counter += 1
    if do_return:
        return counter


def compare_L_deg_pattern(pat_list, seq_list, do_return=False):
    counter = 0
    for i, k in zip(pat_list, seq_list):
        if is_L_deg(i, k):
            counter += 1
    if do_return:
        return counter


def compare_L_bitwize(pat_list, seq_list, do_return=False):
    counter = 0
    for i, k in zip(pat_list, seq_list):
        if bitwize_comp(i,k):
            counter += 1
    if do_return:
        return counter


def assert_proper_comparison_functions(real_letters_list, deg_letters_list_N, deg_letters_list_full):
    simple_comp_res = compare_L_simple_pattern(deg_letters_list_N, real_letters_list, do_return=True)
    simple_comp_bit_res = compare_L_bitwize(deg_letters_list_N, real_letters_list, do_return=True)
    assert(simple_comp_res == simple_comp_bit_res)

    deg_comp_res = compare_L_deg_pattern(deg_letters_list_full, real_letters_list, do_return=True)
    deg_comp_bit_res = compare_L_bitwize(deg_letters_list_full, real_letters_list, do_return=True)
    assert (deg_comp_res == deg_comp_bit_res)


def timing_letters_comparison():
    real_letters_list, deg_letters_list_N, deg_letters_list_full = create_arrays_5_letters_only()
    assert_proper_comparison_functions(real_letters_list, deg_letters_list_N, deg_letters_list_full)

    time_simple_comp = timeit.timeit(lambda: compare_L_simple_pattern(deg_letters_list_N, real_letters_list), number=100000)
    time_simple_comp_bit = timeit.timeit(lambda: compare_L_bitwize(deg_letters_list_N, real_letters_list), number=100000)

    time_deg_comp = timeit.timeit(lambda: compare_L_deg_pattern(deg_letters_list_full, real_letters_list), number=100000)
    time_deg_comp_bit = timeit.timeit(lambda: compare_L_bitwize(deg_letters_list_full, real_letters_list), number=100000)

    print("Compare simple (5 letter) pattern to a sequence timing: ")
    print("Integer comparison takes %.4f seconds" % time_simple_comp)
    print("Bitwise comparison takes %.4f seconds" % time_simple_comp_bit)
    print()
    print("Compare complex (10 degenerate letters) pattern to a sequence timing: ")
    print("Integer comparison takes %.4f seconds" % time_deg_comp)
    print("Bitwise comparison takes %.4f seconds" % time_deg_comp_bit)

    # without numba:
    # in case of simple pattern - integer comparison works ~5 times faster than bitwise
    # in case of complex pattern - integer comparison works ~8 times faster than bitwise
    # with numba:
    # in case of simple pattern - bitwise comparison works ~1.4 times faster than integer
    # in case of complex pattern - integer comparison works ~1.1 times faster than bitwise



def main():
    timing_letters_comparison()


if __name__ == "__main__":
    main()


