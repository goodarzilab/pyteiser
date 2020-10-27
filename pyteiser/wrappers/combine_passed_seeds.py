import argparse
import os
import numpy as np

import IO

def handler(raw_args = None):
    parser = argparse.ArgumentParser()

    parser.add_argument("--passed_seed_folder", help="", type=str)
    parser.add_argument("--passed_seed_filename_template", help="", type=str)
    parser.add_argument("--passed_profiles_folder", help="", type=str)
    parser.add_argument("--passed_profiles_filename_template", help="", type=str)

    parser.add_argument("--combined_seeds_filename", help="output file", type=str)
    parser.add_argument("--combined_profiles_filename", help="output file", type=str)

    parser.add_argument("--input_indices_list_file", help="input: list of indices of files to process", type=str)

    parser.add_argument("--indices_mode", help="compression in the index mode", type=bool)
    parser.add_argument("--index_bit_width", help="number of bits per one index when compressing", type=int)


    parser.set_defaults(
        passed_seed_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_seed/passed_seed_4-7_4-9_4-6_14-20_individual/passed_seed_30k',
        passed_profiles_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_individual/passed_profiles_30k',
        passed_seed_filename_template='passed_seed_4-7_4-9_4-6_14-20_30k',
        passed_profiles_filename_template='passed_profiles_4-7_4-9_4-6_14-20_30k',

        combined_seeds_filename='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_seed/passed_seed_4-7_4-9_4-6_14-20_combined/seeds_passed_100k_tarbp2_utrs.bin',
        combined_profiles_filename='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/profiles_passed_100k_tarbp2_utrs.bin',

        #input_indices_list_file='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/are_profiles_complete/test_2_files.txt',
        input_indices_list_file='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/are_profiles_complete/full_range_2277.txt',

        indices_mode=False,
        index_bit_width = 24,

    )

    args = parser.parse_args(raw_args)

    return args


def get_list_files(passed_seed_folder, passed_seed_filename_template,
                   passed_profiles_folder, passed_profiles_filename_template,
                   input_indices_list_file):
    filenames_list = []

    with open(input_indices_list_file, 'r') as rf:
        full_string = rf.read()
        full_string = full_string.rstrip()
    indices_list_str = full_string.split(', ')
    indices_list = [int(x) for x in indices_list_str]

    for i in sorted(indices_list):
        passed_seed_filename_short = "%s_%d.bin" % (passed_seed_filename_template, i)
        passed_profiles_filename_short = "%s_%d.bin" % (passed_profiles_filename_template, i)
        passed_seed_filename_full = os.path.join(passed_seed_folder, passed_seed_filename_short)
        passed_profiles_filename_full = os.path.join(passed_profiles_folder, passed_profiles_filename_short)
        tuple_filenames = (passed_seed_filename_full, passed_profiles_filename_full)
        filenames_list.append(tuple_filenames)

    return filenames_list


def collect_all_the_passing_seeds(filenames_list, indices_mode):
    all_seeds_passed = []
    all_profiles_passed = []
    for tuple_fn in filenames_list:
        passed_seed_fn, passed_profiles_fn = tuple_fn
        current_seeds_passed = IO.read_seed_pass_individual_file(passed_seed_fn)
        current_profiles_passed = IO.read_profile_pass_individual_file(passed_profiles_fn, indices_mode)
        all_seeds_passed += current_seeds_passed
        if len(current_profiles_passed) > 0: # skip the empty files
            all_profiles_passed.append(current_profiles_passed)

    all_profiles_passed_array = np.concatenate(all_profiles_passed, axis=0)

    return all_seeds_passed, all_profiles_passed_array




def main(raw_args = None):
    args = handler(raw_args)
    filenames_list = get_list_files(args.passed_seed_folder, args.passed_seed_filename_template,
                                    args.passed_profiles_folder, args.passed_profiles_filename_template,
                                    args.input_indices_list_file)
    seeds_passed_list, profiles_passed_array = collect_all_the_passing_seeds(filenames_list, args.indices_mode)
    IO.write_list_of_seeds(seeds_passed_list, args.combined_seeds_filename)
    IO.write_array_of_profiles(profiles_passed_array, args.combined_profiles_filename,
                               args.indices_mode, args.index_bit_width)


if __name__ == "__main__":
    main()



