import argparse


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--incomplete_files_summary", help="input: list of incomplete profile files in an unparsed format", type=str)
    parser.add_argument("--files_indices", help="output: list of indices of files to process", type=str)
    parser.set_defaults(
        incomplete_files_summary='/Users/student/Documents/hani/programs/pyteiser/data/testing_data/are_profiles_complete/profiles_tarbp2.txt',
        files_indices='/Users/student/Documents/hani/programs/pyteiser/data/testing_data/are_profiles_complete/indices_incomplete_tarbp2.txt',
    )

    args = parser.parse_args()

    return args


def parse_list_incomplete_files(input_file, output_file):
    file_numbers_list = []
    with open(input_file, 'r') as rf:
        for i, line in enumerate(rf):
            if not line.startswith("File "):
                continue
            long_filename = line.split(' ')[1]
            short_filename = long_filename.split('/')[-1]
            file_number_string = short_filename.split('_')[-1].replace('.bin','')
            file_number = int(file_number_string)
            file_numbers_list.append(file_number)

    strings_file_numbers_list = [str(x) for x in sorted(file_numbers_list)]
    string_to_write = ", ".join(strings_file_numbers_list)

    with open(output_file, 'w') as wf:
        wf.write(string_to_write)


def main():
    args = handler()
    parse_list_incomplete_files(args.incomplete_files_summary,
                                                    args.files_indices)


if __name__ == "__main__":
    main()

