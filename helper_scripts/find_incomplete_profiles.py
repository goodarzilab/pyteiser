import argparse
import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
sys.path.append("/wynton/home/goodarzi/khorms/pyteiser_root/pyteiser/pyteiser")

import IO


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--profiles_folder", help="", type=str)
    parser.add_argument("--outfile", help="", type=str)
    parser.add_argument("--expected_length", help="", type=int)
    parser.add_argument("--indices_mode", help="compression in the index mode", type=bool)
    parser.set_defaults(
        #profiles_folder = '/wynton/home/goodarzi/khorms/pyteiser_root/data/profiles/profiles_4-7_4-9_4-6_14-20/profiles_per_file_30k',
        profiles_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/profiles/SNRNPA1/profiles_per_file_30k',
        #outfile = '/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/are_profiles_complete/profiles_tarbp2.txt',
        outfile='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/are_profiles_complete/profiles_snrnp1.txt',
        expected_length = 30000,
        indices_mode = False,
        identify_broken_ones_only = False,
    )

    args = parser.parse_args()

    return args


def get_list_profile_files(profiles_folder):
    short_filenames_list = os.listdir(profiles_folder)
    short_filenames_list = [x for x in short_filenames_list if x.endswith('.bin')]
    long_filenames_list = [os.path.join(profiles_folder, x) for x in short_filenames_list]
    return sorted(long_filenames_list)



def report_incomplete_profiles(list_profile_filenames, output_file, expected_length,
                               indices_mode = False, identify_broken_ones_only = False):
    any_files_incomplete = False
    with open(output_file, 'w') as wf:
        counter = 0
        for fn in list_profile_filenames:
            with open(fn, 'rb') as rf:
                bitstring = rf.read()
                is_broken = False
                if indices_mode:
                    try:
                        decompressed_profiles_array = IO.decompress_profiles_indices(bitstring)
                    except:
                        decompressed_profiles_array = []
                        is_broken = True
                else:
                    try:
                        decompressed_profiles_array = IO.decompress_profiles(bitstring)
                    except:
                        decompressed_profiles_array = []
                        is_broken = True
                if not identify_broken_ones_only:
                    if len(decompressed_profiles_array) != expected_length:
                        any_files_incomplete = True
                        string_to_write = 'File %s contains %d profiles instead of %d\n' % \
                                        (fn, len(decompressed_profiles_array), expected_length)
                        wf.write(string_to_write)
                else:
                    if is_broken:
                        string_to_write = 'File %s contains is broken\n' % fn
                        wf.write(string_to_write)


            counter += 1
            # print("Processed file number ", counter)

        if any_files_incomplete:
            wf.write("Some incomplete files have been identified!\n")
        wf.write("All the files in the specified folder have been processed")


def main():
    args = handler()
    list_profile_filenames = get_list_profile_files(args.profiles_folder)
    report_incomplete_profiles(list_profile_filenames, args.outfile, args.expected_length,
                               args.indices_mode, args.identify_broken_ones_only)


if __name__ == "__main__":
    main()

