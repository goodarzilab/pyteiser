import argparse
import os
import sys
import numpy as np


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--passed_seed_folder", help="", type=str)
    parser.add_argument("--passed_seed_filename_template", help="", type=str)
    parser.add_argument("--passed_profiles_folder", help="", type=str)
    parser.add_argument("--passed_profiles_filename_template", help="", type=str)

    parser.add_argument("--combined_seeds_filename", help="output file", type=str)
    parser.add_argument("--combined_profiles_filename", help="output file", type=str)

    parser.add_argument("--input_indices_list_file", help="input: list of indices of files to process", type=str)


    parser.set_defaults(
        passed_seed_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_seed/passed_seed_4-7_4-9_4-6_14-20_individual/passed_seed_30k',
        passed_profiles_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_individual/passed_profiles_30k',
        passed_seed_filename_template='passed_seed_4-7_4-9_4-6_14-20_30k',
        passed_profiles_filename_template='passed_profiles_4-7_4-9_4-6_14-20_30k',

        combined_seeds_filename='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_seed/passed_seed_4-7_4-9_4-6_14-20_combined/seeds_passed_100k_tarbp2_utrs.bin',
        combined_profiles_filename='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/profiles_passed_100k_tarbp2_utrs.bin',

        input_indices_list_file='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/are_profiles_complete/test_2_files.txt',
        #input_indices_list_file='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/are_profiles_complete/full_range_2277.txt',

    )

    args = parser.parse_args()

    return args


def import_modules():
    current_script_path = sys.argv[0]
    package_home_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if package_home_path not in sys.path:
        sys.path.append(package_home_path)

    global structures
    global IO

    import structures
    import IO


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


def collect_all_the_passing_seeds(filenames_list):
    all_seeds_passed = []
    all_profiles_passed = []
    for tuple_fn in filenames_list:
        passed_seed_fn, passed_profiles_fn = tuple_fn
        current_seeds_passed = IO.read_seed_pass_individual_file(passed_seed_fn)
        current_profiles_passed = IO.read_profile_pass_individual_file(passed_profiles_fn)
        all_seeds_passed += current_seeds_passed
        all_profiles_passed.append(current_profiles_passed)

    all_profiles_passed_array = np.concatenate(all_profiles_passed, axis=0)

    return all_seeds_passed, all_profiles_passed_array


def write_seeds_passed(seeds_passed_list, combined_seeds_filename):
    seeds_bitstrings = []

    for motif in seeds_passed_list:
        motif.compress()
        seeds_bitstrings.append(motif.bytestring)

    total_bitstring = b''.join(seeds_bitstrings)
    #print("Total number is: ", len(seeds_bitstrings))

    with open(combined_seeds_filename, 'wb') as wf:
        wf.write(total_bitstring)


def write_profiles_passed(profiles_passed_array, combined_profiles_filename):
    with open(combined_profiles_filename, 'wb') as wf:
        for i in range(profiles_passed_array.shape[0]):
            current_profile = structures.w_profile(profiles_passed_array[i].shape[0])
            current_profile.values = profiles_passed_array[i]
            current_profile.compress()
            wf.write(current_profile.bytestring)


def main():
    import_modules()
    args = handler()
    filenames_list = get_list_files(args.passed_seed_folder, args.passed_seed_filename_template,
                                    args.passed_profiles_folder, args.passed_profiles_filename_template,
                                    args.input_indices_list_file)
    seeds_passed_list, profiles_passed_array = collect_all_the_passing_seeds(filenames_list)
    write_seeds_passed(seeds_passed_list, args.combined_seeds_filename)
    write_profiles_passed(profiles_passed_array, args.combined_profiles_filename)


if __name__ == "__main__":
    main()



