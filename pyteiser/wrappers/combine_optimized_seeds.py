import numpy as np
import argparse
import math

import os
import sys




def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task_mapping_file", help="", type=str)

    parser.add_argument("--optimized_seeds_filename", help="output: optimized seeds", type=str)
    parser.add_argument("--optimized_profiles_filename", help="output: profiles of optimized seeds", type=str)
    parser.add_argument("--optimized_MI_pv_zscores_filename", help="output: MI values, p-values and z-scores "
                                                                   "of optimized seeds", type=str)
    parser.add_argument("--robustness_array_filename", help="output: vector indicating which seeds "
                                                            "have passed the robustness test", type=str)

    parser.add_argument("--print_qstat", help="", type=str)
    parser.add_argument("--path_to_qstat", help="", type=str)


    parser.set_defaults(
        optimized_seeds_filename='/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds_optimized.bin',
        optimized_profiles_filename='/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_optimized.bin',
        optimized_MI_pv_zscores_filename='/Users/student/Documents/hani/programs/pyteiser/data/optimized_seeds_characteristics/seeds_4-7_4-9_4-6_14-20_individual/test_1_2_characteristics.bin',
        robustness_array_filename='/Users/student/Documents/hani/programs/pyteiser/data/seeds_robustness/seeds_4-7_4-9_4-6_14-20_individual/test_1_2_robustness.bin',

        path_to_qstat='/opt/sge/bin/lx-amd64/qstat',
        print_qstat='y',

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



















