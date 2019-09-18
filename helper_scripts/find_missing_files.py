import os
import argparse


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--folder", help="folder where to look for missing", type=str)
    parser.set_defaults(
    )

    args = parser.parse_args()

    return args


def find_all_indices(inp_folder):
    all_filenames = os.listdir(inp_folder)
    all_filenames = [x for x in all_filenames if '.bin' in x]
    all_indices_string = [x.split('_')[-1].replace('.bin', '') for x in all_filenames]
    all_indices_int = [int(x) for x in all_indices_string]
    return all_indices_int


def find_missing_indices(all_indices_list, min_index = 1):
    max_index = max(all_indices_list)
    print("Maximal index found is ", max_index)
    full_set_indices = set(list(range(min_index, max_index)))
    existing_set_indices = set(all_indices_list)
    difference = full_set_indices.difference(existing_set_indices)
    missing_ones_list = sorted([int(x) for x in list(difference)])
    missing_ones_string = ", ".join(missing_ones_list)
    print("The missing indices are: ", missing_ones_string)


def main():
    args = handler()
    all_indices_list = find_all_indices(args.folder)
    find_missing_indices(all_indices_list)


if __name__ == "__main__":
    main()