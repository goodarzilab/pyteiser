import argparse
import os
import sys
import numpy as np


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--combined_seeds_filename", help="output file", type=str)
    parser.add_argument("--combined_profiles_filename", help="output file", type=str)

    parser.add_argument("--exp_mask_file", help="file with binary expression file, pre-overlapped with "
                                                "the reference transcriptome", type=str)
    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)

    parser.set_defaults(
        # combined_seeds_filename='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_seed/passed_seed_4-7_4-9_4-6_14-20_combined/seeds_passed_100k_tarbp2_utrs.bin',
        # combined_profiles_filename='/wynton/home/goodarzi/khorms/pyteiser_root/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/profiles_passed_100k_tarbp2_utrs.bin',
        # exp_mask_file='/wynton/home/goodarzi/khorms/pyteiser_root/data/mask_files/TARBP2_decay_t_score_mask.bin',
        combined_seeds_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_seeds/passed_seed_4-7_4-9_4-6_14-20_combined/test_1_2_seeds.bin",
        combined_profiles_filename="/Users/student/Documents/hani/programs/pyteiser/data/passed_profiles/passed_profiles_4-7_4-9_4-6_14-20_combined/test_1_2_profiles.bin",
        exp_mask_file='/Users/student/Documents/hani/programs/pyteiser/data/mask_files/TARBP2_decay_t_score_mask.bin',

        nbins=15,

    )

    args = parser.parse_args()

    return args


def import_modules():
    current_script_path = sys.argv[0]
    package_home_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if package_home_path not in sys.path:
        sys.path.append(package_home_path)

    global MI
    global IO

    import MI
    import IO



def min_cond_info_normalized():
    for i in range(len(accepted_seeds_list)):
        ith_accepted_profile_full = profiles_passed[accepted_seeds_list[i]]
        ith_accepted_profile = ith_accepted_profile_full[index_array]
        combined_profile = np.logical_or(curr_accepted_profile, profile_being_analyzed)

        MI_accepted_vs_expr = MI.mut_info(ith_accepted_profile, discr_exp_profile, x_bins=2, y_bins=2)
        MI_combined_vs_expr = MI.mut_info(combined_profile, discr_exp_profile, x_bins=2, y_bins=2)
        cond_MI = MI_combined_vs_expr - MI_accepted_vs_expr





        # calculate conditional information I(Cnew;E|Cexist)  = I(A;E|B)  =  I(A,B;E) - I(A;E)
        M_q = n-th profile among the ones we've kept
        B_q = the profile we are analyzing
        AB_q = combineQuantizedVectors(M_q, B_q, n, DA, DB); - pairwise OR of M_q and B_q
        mi_ab_e = CalculateMIbasic(AB_q, E_q, n, DAB, DE);
        mi_a_e = CalculateMIbasic(M_q, E_q, n, DA, DE);
        cmi = mi_ab_e - mi_a_e;

        # calculate MI I(Cnew;Cexist)
        mi_a_b = CalculateMIbasic(M_q, B_q, n, DA, DB);

    return



def filter_CMI(seeds_passed, profiles_passed,
               discr_exp_profile, index_array):
    indices_accepted_profiles = []
    for i in range(len(seeds_passed)):
        seed = seeds_passed[i]
        profile_full = profiles_passed[i]
        profile_being_analyzed = profile_full[index_array]
        minratio = minCondInfoNormalized(opt_motifs, mbins, accepted_seeds_count, M_q, mbins, E_q, ebins, seq_count,
                                         minr, & midx, sequences, t_seq_count, h_rna_ind, dG_t)

        min_cond_info_normalized()

        # current_MI = MI.mut_info(active_profile, discr_exp_profile_tarbp2, x_bins=2, y_bins=15)
        # pvalue, z_score = statistic_tests.MI_get_pvalue_and_zscore(active_profile, discr_exp_profile_tarbp2,
        #                                                            15, current_MI, n_permutations=1000)



def main():
    import_modules()
    args = handler()

    index_array, values_array = IO.unpack_mask_file(args.exp_mask_file)
    discr_exp_profile = MI.discretize_exp_profile(index_array, values_array, nbins = args.nbins)
    seeds_passed = IO.read_motif_file(args.combined_seeds_filename)
    profiles_passed = IO.unpack_profiles_file(args.combined_profiles_filename)

    print(profiles_passed.shape)



if __name__ == "__main__":
    main()