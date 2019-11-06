import os
import sys
import numpy as np
import argparse
import math
import subprocess


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--folded_bin_filename", type=str)


    parser.set_defaults(
        folded_bin_filename = '/Users/student/Documents/hani/programs/pyteiser/data/transcriptome_folding/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_folded.bin',
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
    global glob_var

    import structures
    import IO
    import type_conversions
    import glob_var



def main():
    import_modules()
    args = handler()

    arrays_list = IO.read_multiple_np_arrays(args.folded_bin_filename,
                                             dtype = np.dtype('uint8'))
    for array in arrays_list:
        print(array.shape)
        print(array)


if __name__ == "__main__":
    main()
