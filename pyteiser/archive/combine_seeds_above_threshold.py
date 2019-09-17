import argparse
import os
import sys
import numpy as np

current_script_path = sys.argv[0]
sys.path.append("/wynton/home/goodarzi/khorms/pyteiser_root/pyteiser/pyteiser")

import IO


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--thresholds_folder", help="", type=str)
    parser.add_argument("--seed_folder", help="", type=str)
    parser.add_argument("--MI_values_folder", help="", type=str)
    parser.add_argument("--thresholds_filename_template", help="", type=str)
    parser.add_argument("--seed_filename_template", help="", type=str)

    parser.add_argument("--input_indices_list_file", help="input: list of indices of files to process", type=str)
    parser.add_argument("--seeds_passed_filename", help="output file", type=str)

    parser.set_defaults(
        thresholds_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/thresholds/thresholds_4-7_4-9_4-6_14-20/thresholds_per_file_30k/permutations_100k',
        seed_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/seeds/seeds_4-7_4-9_4-6_14-20/motifs_per_file_30k',
        MI_values_folder='/wynton/home/goodarzi/khorms/pyteiser_root/data/MI_values/MI_values_4-7_4-9_4-6_14-20/MI_values_per_file_30k',
        thresholds_filename_template='thresholds_4-7_4-9_4-6_14-20_30k',
        MI_values_filename_template='MI_values_4-7_4-9_4-6_14-20_30k',
        seed_filename_template='seeds_4-7_4-9_4-6_14-20_30k',

        #input_indices_list_file='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/are_profiles_complete/test_2_files.txt',
        input_indices_list_file='/wynton/home/goodarzi/khorms/pyteiser_root/testing_data/are_profiles_complete/full_range_2277.txt',
        seeds_passed_filename='/wynton/home/goodarzi/khorms/pyteiser_root/data/thresholds/thresholds_4-7_4-9_4-6_14-20/thresholds_per_file_30k/seeds_past_threshold_p100k/seeds_passed_100k_tarbp2_utrs.bin',

    )

    args = parser.parse_args()

    return args


def get_list_files(thresholds_folder, thresholds_filename_template,
                   seed_folder, seed_filename_template,
                   MI_values_folder, MI_values_filename_template,
                   input_indices_list_file):
    filenames_tuples_list = []

    with open(input_indices_list_file, 'r') as rf:
        full_string = rf.read()
        full_string = full_string.rstrip()
    indices_list_str = full_string.split(', ')
    indices_list = [int(x) for x in indices_list_str]

    for i in sorted(indices_list):
        thresholds_filename_short = "%s_%d.bin" % (thresholds_filename_template, i)
        thresholds_filename_full = os.path.join(thresholds_folder, thresholds_filename_short)
        seed_filename_short = "%s_%d.bin" % (seed_filename_template, i)
        seed_filename_full = os.path.join(seed_folder, seed_filename_short)
        MI_values_filename_short = "%s_%d.bin" % (MI_values_filename_template, i)
        MI_values_filename_full = os.path.join(MI_values_folder, MI_values_filename_short)

        filenames_tuple = (thresholds_filename_full, seed_filename_full, MI_values_filename_full)
        filenames_tuples_list.append(filenames_tuple)

    return filenames_tuples_list


def extract_seeds_by_threshold(thresh_fn, seed_fn, mi_values_fn):
    threshold_value = IO.read_seed_significancy_threshold(thresh_fn)
    w_motifs_list = IO.read_motif_file(seed_fn)
    MI_values_array, nbins = IO.read_MI_values(mi_values_fn)

    seed_indices_sorted = np.argsort(MI_values_array)[::-1]
    indices_passed = seed_indices_sorted[0 : threshold_value]
    w_motifs_passed = [w_motifs_list[x] for x in indices_passed]

    return w_motifs_passed



def collect_all_the_passing_seeds(filenames_pairs_list):
    all_seeds_passed = []
    for tup in filenames_pairs_list:
        thresh_fn, seed_fn, mi_values_fn = tup
        current_seeds_passed = extract_seeds_by_threshold(thresh_fn, seed_fn, mi_values_fn)
        all_seeds_passed += current_seeds_passed

    return all_seeds_passed


def write_seeds_passed(seeds_passed_list, seeds_passed_filename):
    seeds_bitstrings = []

    for motif in seeds_passed_list:
        motif.compress()
        seeds_bitstrings.append(motif.bytestring)

    total_bitstring = b''.join(seeds_bitstrings)

    with open(seeds_passed_filename, 'wb') as wf:
        wf.write(total_bitstring)




def main():
    args = handler()
    filenames_pairs_list = get_list_files(args.thresholds_folder, args.thresholds_filename_template,
                                          args.seed_folder, args.seed_filename_template,
                                          args.MI_values_folder, args.MI_values_filename_template,
                                          args.input_indices_list_file)
    seeds_passed_list = collect_all_the_passing_seeds(filenames_pairs_list)
    write_seeds_passed(seeds_passed_list, args.seeds_passed_filename)



if __name__ == "__main__":
    main()



