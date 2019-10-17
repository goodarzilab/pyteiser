import numpy as np
import argparse
import math

import os
import sys


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--combined_seeds_filename", help="", type=str)
    parser.add_argument("--combined_profiles_filename", help="", type=str)
    parser.add_argument("--combined_MI_pv_zscores_filename", help="", type=str)
    parser.add_argument("--combined_robustness_filename", help="", type=str)

    parser.set_defaults(
        combined_seeds_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/tarbp2/seed_optimized_100k_tarbp2_utrs_10k.bin',
        combined_profiles_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/tarbp2/profiles_optimized_100k_tarbp2_utrs_10k.bin',
        combined_MI_pv_zscores_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/tarbp2/seed_characteristics_optimized_100k_tarbp2_utrs_10k.bin',
        combined_robustness_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/tarbp2/robustness_optimized_100k_tarbp2_utrs_10k.bin',

        # combined_seeds_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/snrnpa1/seed_optimized_100k_snrnpa1_10k.bin',
        # combined_profiles_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/snrnpa1/profiles_optimized_100k_snrnpa1_10k.bin',
        # combined_MI_pv_zscores_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/snrnpa1/seed_characteristics_optimized_100k_snrnpa1_10k.bin',
        # combined_robustness_filename='/Users/student/Documents/hani/programs/pyteiser/data/combined_optimized_seeds/snrnpa1/robustness_optimized_100k_snrnpa1_10k.bin',
    )

    args = parser.parse_args()

    return args



def import_modules():
    package_home_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if package_home_path not in sys.path:
        sys.path.append(package_home_path)

    global structures
    global IO
    global type_conversions

    import structures
    import IO



def read_seeds_and_characteristics(args):
    seeds_optimized = IO.read_motif_file(args.combined_seeds_filename)
    profiles_optimized = IO.unpack_profiles_file(args.combined_profiles_filename)
    seed_charact_array = IO.read_np_array(args.combined_MI_pv_zscores_filename, np.dtype('float64'))
    robustness_array = IO.read_np_array(args.combined_robustness_filename, np.dtype('bool'))

    return seeds_optimized, profiles_optimized, \
           seed_charact_array, robustness_array


def print_read_out_arrays(seeds_optimized, profiles_optimized,
                          seed_charact_array, robustness_array):
    MI_values_array = seed_charact_array[:, 0]
    seed_indices_sorted = np.argsort(MI_values_array)[::-1]

    for i in seed_indices_sorted:
        is_robust = robustness_array[i]
        if not is_robust:
            continue

        string_to_print = "Seed %d. Sequence: %s. " % (i, seeds_optimized[i].print_sequence(return_string=True))
        string_to_print += "MI: %.3f; p-value: %.5f; z-score: %.1f. " % \
                           (seed_charact_array[i,0], seed_charact_array[i,1], seed_charact_array[i,2])
        string_to_print += "It binds %d sequences." % (profiles_optimized[i].sum())

        print(string_to_print)




def main():
    import_modules()
    args = handler()

    seeds_optimized, profiles_optimized, \
    seed_charact_array, robustness_array = read_seeds_and_characteristics(args)


    print_read_out_arrays(seeds_optimized, profiles_optimized,
                          seed_charact_array, robustness_array)





if __name__ == "__main__":
    main()