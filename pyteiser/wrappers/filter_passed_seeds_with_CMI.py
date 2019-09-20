import argparse
import os
import sys
import numpy as np


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--combined_seeds_filename", help="file with seeds that passed the thresholding phase", type=str)
    parser.add_argument("--combined_profiles_filename", help="file with profiles of seeds "
                                                             "that passed the thresholding phase", type=str)

    parser.add_argument("--unique_seeds_filename", help="output: best representatives of each family", type=str)
    parser.add_argument("--unique_profiles_filename", help="output: profiles of best representatives of each family",
                                                                                                        type=str)
    parser.add_argument("--families_classification_filename", help="output: classification of all the passed seeds"
                                                                   "to unique families", type=str)

    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)
    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)
    parser.add_argument("--min_ratio", help="threshold on ratio of CMI/MI for the conditional "
                                            "information test for seed novelty", type=int)
    parser.add_argument("--do_print", help="number of bins for discretization of expression profile", type=int)

    parser.set_defaults(
        # combined_seeds_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/seeds_passed_100k_tarbp2_utrs.bin",
        # combined_profiles_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/profiles_passed_100k_tarbp2_utrs.bin",

        combined_seeds_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds.bin",
        combined_profiles_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles.bin",

        unique_seeds_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds_unique.bin",
        unique_profiles_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles_unique.bin",
        families_classification_filename="/Users/student/Documents/hani/programs/pyteiser/data/seeds_family_classification/seeds_4-7_4-9_4-6_14-20_combined/test_1_2_classification.bin",

        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',



        nbins=15,
        min_ratio=5,
        do_print=False,

    )

    args = parser.parse_args()

    return args


def import_modules():
    package_home_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if package_home_path not in sys.path:
        sys.path.append(package_home_path)

    global MI
    global IO

    import MI
    import IO



# loop through existing profiles. For each existing profile:
# calculate conditional information I(new profile ; expression | existing profile)
# calculate I(new profile ; existing profile)
# take their ratio
# if a new profile is almost the same as an old profile, the MI between them two will be large
# and the CMI between the new profile and the expression profile given existing profile will be low
# so we want to catch the cases like that and filter them out
def min_CI_normalized_test(counter, accepted_seeds_list, profiles_passed,
                             discr_exp_profile, nbins, index_array, min_ratio,
                           do_print=False):
    profile_full = profiles_passed[counter]
    profile_being_analyzed = profile_full[index_array]

    for i in range(len(accepted_seeds_list)):
        ith_accepted_profile_full = profiles_passed[accepted_seeds_list[i]]
        ith_accepted_profile = ith_accepted_profile_full[index_array]

        cond_inf = MI.cond_mut_info(profile_being_analyzed, discr_exp_profile, ith_accepted_profile,
                      x_bins=2, y_bins=nbins, z_bins=2)
        mut_inf = MI.mut_info(profile_being_analyzed, ith_accepted_profile, x_bins=2, y_bins=2)
        if np.isclose(mut_inf, 0., atol=1e-16):
            mut_inf = 1e-16
        ratio = cond_inf / mut_inf

        print("Comparing seed #%d to an existing seed #%d. The ratio is %.2f" % (counter, i, ratio))

        if ratio < min_ratio:
            return False, i # return index of accepted seed that is similar to the current one
    return True, 0



def filter_CMI(seeds_passed, profiles_passed,
               discr_exp_profile, index_array,
               nbins, min_ratio, do_print = False):
    classification_array = np.zeros(len(seeds_passed), dtype=np.int32)
    indices_accepted_profiles = [0]
    for k in range(1,len(seeds_passed)):
        # seed = seeds_passed[i]
        is_seed_new, similar_one_index = min_CI_normalized_test(k, indices_accepted_profiles,
                                profiles_passed, discr_exp_profile, nbins, index_array,
                                min_ratio, do_print = do_print)
        if is_seed_new:
            classification_array[k] = len(indices_accepted_profiles)
            indices_accepted_profiles.append(k)
            if do_print:
                print(indices_accepted_profiles)
                print('\n')
        else:
            classification_array[k] = similar_one_index

    accepted_profiles_array = np.array(indices_accepted_profiles)
    seeds_unique = [seeds_passed[x] for x in accepted_profiles_array]
    profiles_unique = profiles_passed[accepted_profiles_array]

    return seeds_unique, profiles_unique


def calculate_MIs_all_seeds():
    active_profile = profile[index_array]
    MI_values_array[i] = MI.mut_info(active_profile, discr_exp_profile, x_bins=2, y_bins=nbins)



def main():
    import_modules()
    args = handler()

    index_array, values_array = IO.unpack_mask_file(args.exp_mask_file)
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins = args.nbins)
    seeds_passed = IO.read_motif_file(args.combined_seeds_filename)
    profiles_passed = IO.unpack_profiles_file(args.combined_profiles_filename)

    seeds_unique, profiles_unique = filter_CMI(seeds_passed, profiles_passed,
                                               discr_exp_profile, index_array,
                                               args.nbins, args.min_ratio,
                                               do_print = True)
    # IO.write_list_of_seeds(seeds_unique, args.unique_seeds_filename)
    # IO.write_array_of_profiles(profiles_unique, args.unique_profiles_filename)



if __name__ == "__main__":
    main()