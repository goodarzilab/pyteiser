import numpy as np
import argparse

import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import IO
import matchmaker
import type_conversions


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--profiles", type=str)


    parser.set_defaults(
        profiles="/Users/student/Documents/hani/programs/pyteiser/data/test_profiles/test_motifs_101.bin"
    )

    args = parser.parse_args()

    return args




def main():
    args = handler()


if __name__ == "__main__":
    main()
