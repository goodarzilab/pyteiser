import argparse
import os
import sys

current_script_path = sys.argv[0]
sys.path.append("/wynton/home/goodarzi/khorms/pyteiser_root/pyteiser/pyteiser")

import IO


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--thresholds_folder", help="", type=str)
    parser.add_argument("--outfile", help="", type=str)
    parser.set_defaults(
        thresholds_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/thresholds/thresholds_4-7_4-9_4-6_14-20/thresholds_per_file_30k/permutations_100k',
        outfile='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/threshold_values/Sept_10_tarbp_thresh_values.txt',

    )

    args = parser.parse_args()

    return args


def get_list_files(profiles_folder):
    short_filenames_list = os.listdir(profiles_folder)
    short_filenames_list = [x for x in short_filenames_list if x.endswith('.bin')]
    long_filenames_list = [os.path.join(profiles_folder, x) for x in short_filenames_list]
    return sorted(long_filenames_list)



def write_threshold_values(list_profile_filenames, output_file):
    with open(output_file, 'w') as wf:
        for fn in list_profile_filenames:
            current_threshold_value = IO.read_seed_significancy_threshold(fn)
            wf.write("{}\t{}\n".format(os.path.basename(fn), current_threshold_value))


def main():
    args = handler()
    list_profile_filenames = get_list_files(args.thresholds_folder)
    write_threshold_values(list_profile_filenames, args.outfile)


if __name__ == "__main__":
    main()


