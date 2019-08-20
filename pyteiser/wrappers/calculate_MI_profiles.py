import numpy as np
import argparse

import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import IO
import MI
import matchmaker
import type_conversions


MASK_OUT_SEED_VALUE = np.float64(-1)


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--profiles_bin_file", help="file with occurence profiles", type=str)
    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)
    parser.add_argument("--MI_values_file", help="output file where calculated MI values are stored", type=str)

    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)
    parser.add_argument("--min_occurences", help="minimal number of seed occurence in the transcriptome"
                                                 " for a seed to be considered", type=int)


    parser.set_defaults(
        profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/test_motifs_101.bin",
        #profiles_bin_file="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/profiles_4-7_4-9_4-6_14-20_30k_1.bin",
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',
        MI_values_file='/Users/student/Documents/hani/programs/pyteiser/data/MI_values/MI_test_motifs_101.bin',

        nbins = 3,
        min_occurences = 5,
    )

    args = parser.parse_args()

    return args


def unpack_profiles_and_mask(args, do_print=False):
    with open(args.profiles_bin_file, 'rb') as rf:
        bitstring = rf.read()
        decompressed_profiles_array = IO.decompress_profiles(bitstring)
        if do_print:
            print("%d profiles have been loaded" % len(decompressed_profiles_array))

    with open(args.exp_mask_file, 'rb') as rf:
        bitstring = rf.read()
        index_array, values_array = IO.decompress_exp_mask_file(bitstring)
        if do_print:
            print("Expression values are provided for %d out of %d transcripts in the reference transcriptome" %
                  (index_array.sum(), index_array.shape[0]))

    try:
        assert (decompressed_profiles_array[0].shape[0] == index_array.shape[0])
    except AssertionError:
        print("Error: occurence profiles were calculated for some other reference transcriptome. The length of the "
              "profiles is %d and the length of the transcriptome provided is %d" %
              (decompressed_profiles_array[0].shape[0], index_array.shape[0]))
        sys.exit(1)

    return decompressed_profiles_array, index_array, values_array


def discretize_exp_profile(index_array, values_array, nbins):
    active_values_array = values_array[index_array]
    quant_values_array = MI.discretize(active_values_array, bins=nbins, disc = "equalfreq")
    return quant_values_array



def calculate_MI_for_seeds(decompressed_profiles_array, index_array, discr_exp_profile,
                       min_occurences, do_print=False):
    MI_values_array = np.zeros(index_array.shape[0], dtype=np.float64)

    for i, profile in enumerate(decompressed_profiles_array):
        active_profile = profile[index_array]

        if active_profile.sum() <= min_occurences:
            MI_values_array[i] = MASK_OUT_SEED_VALUE
            continue

        MI_values_array[i] = MI.mut_info(active_profile, discr_exp_profile)

        if do_print:
            if i % 1000 == 0:
                print("Profile number %d has been calculated" % i)

    return MI_values_array


def write_MI_values(MI_values_array, MI_values_file):
    length_uint32 = np.array([MI_values_array.shape[0]], dtype=np.uint32)
    length_bitstring = length_uint32.tobytes()
    values_bytes = MI_values_array.tobytes()
    MI_values_bytestring = length_bitstring + values_bytes

    with open(MI_values_file, 'wb') as wf:
        wf.write(MI_values_bytestring)






def main():
    args = handler()
    decompressed_profiles_array, index_array, values_array = unpack_profiles_and_mask(args, do_print=True)
    discr_exp_profile = discretize_exp_profile(index_array, values_array, args.nbins)
    MI_values_array = calculate_MI_for_seeds(decompressed_profiles_array, index_array, discr_exp_profile,
                                         args.min_occurences, do_print = True)
    write_MI_values(MI_values_array, args.MI_values_file)


    # proceed with line 179 in mi_find_seed.c



if __name__ == "__main__":
    main()
