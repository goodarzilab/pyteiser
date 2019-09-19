import argparse
import os
import sys
import numpy as np


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--combined_seeds_filename", help="output file", type=str)
    parser.add_argument("--combined_profiles_filename", help="output file", type=str)

    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)
    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)

    parser.set_defaults(
        # combined_seeds_filename='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_seed/passed_seed_4-7_4-9_4-6_14-20_combined/seeds_passed_100k_tarbp2_utrs.bin',
        # combined_profiles_filename='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/profiles_passed_100k_tarbp2_utrs.bin',
        # exp_mask_file='/wynton/home/goodarzi/khorms/pyteiser_root/data/mask_files/TARBP2_decay_t_score_mask.bin',
        combined_seeds_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds.bin",
        combined_profiles_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles.bin",
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',

        nbins=15,

    )

    args = parser.parse_args()

    return args






def import_modules():
    current_script_path = sys.argv[0]
    package_home_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if package_home_path not in sys.path:
        sys.path.append(package_home_path)

    global MI
    global IO

    import MI
    import IO



def main():
    import_modules()
    args = handler()

    index_array, values_array = IO.unpack_mask_file(args.exp_mask_file)
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins = args.nbins)
    seeds_passed = IO.read_motif_file(args.combined_seeds_filename)
    profiles_passed = IO.unpack_profiles_file(args.combined_profiles_filename)

    print(profiles_passed.shape)



if __name__ == "__main__":
    main()