import os
import argparse

from .. import IO
from . import binarize_sequences
from . import preprocess_custom_expression_profile
from . import calculate_seed_profiles
from . import calculate_MI_profiles
from . import choose_significant_seeds_v3
from . import filter_passed_seeds_with_CMI
from . import combine_passed_seeds
from . import optimize_seeds_single_chunk
from . import report_csvs


def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rna_fastafile", help="fasta file with RNA sequences", type=str)
    parser.add_argument("--exp_values_file", help="expression values in a csv format", type=str)
    parser.add_argument("--anno_name_column", help="column name in exp_values file that contains annotations", type=str)
    parser.add_argument("--measur_column", help="column name in exp_values file that contains expression measurements", type=str)
    parser.add_argument("--seeds_file", help="file with seeds in binary format", type=str)
    parser.add_argument("--temp_folder", help="folder to write temporary files to", type=str)
    parser.add_argument("--out", help="output folder", type=str)

    parser.add_argument("--nbins", help="number of bins for discretization of expression profile", type=int)
    parser.add_argument("--min_occurences", help="minimal number of seed occurence in the transcriptome"
                                                 " for a seed to be considered", type=int)
    parser.add_argument("--maxfreq", help="maximal seed frequency in the sequences analyzed", type=float)
    parser.add_argument("--n_permutations", help="number of permutations for the rank test for a seed", type=int)
    parser.add_argument("--max_pvalue", help="p-value threshold", type=float)
    parser.add_argument("--min_zscore", help="z-score threshold", type=float)
    parser.add_argument("--step_1_jump", help="step size at the 1st round of empirical threshold search", type=int)
    parser.add_argument("--step_2_min_interval", help="resolution of empirical threshold search", type=int)
    parser.add_argument("--step_1_min_fraction", help="minimal fraction of passing seeds for the 1st round of empirical threshold search", type=float)
    parser.add_argument("--step_2_min_fraction", help="minimal fraction of passing seeds for the 2nd round of empirical threshold search", type=float)
    parser.add_argument("--step_3_min_fraction", help="minimal fraction of passing seeds for the 3rd round of empirical threshold search", type=float)
    parser.add_argument("--min_ratio", help="threshold on ratio of CMI/MI for the conditional "
                                            "information test for seed novelty", type=int)

    parser.add_argument("--are_input_seeds_degenerate", help="do input seeds contain degenerate nucleotides", type=bool)
    parser.add_argument("--indices_mode", help="compression in the index mode", type=bool)
    parser.add_argument("--index_bit_width", help="number of bits per one index when compressing", type=int)
    parser.add_argument("--random_noseed", help="when choosing the order of positions to optimize, "
                                                "do not set the random number generator to a specific seed", type=bool)
    parser.add_argument("--jackknife_n_permutations", help="number of permutations for pvalue calculation in "
                                                           "jackknife test", type=int)
    parser.add_argument("--jackknife_max_pvalue", help="maximal pvalue for jackknife test", type=float)
    parser.add_argument("--jackknife_n_samples", help="how many permutations to do in jackknife test", type=int)
    parser.add_argument("--jackknife_fraction_retain", help="what fraction of the sample to retain for each test",
                                                                                    type=float)
    parser.add_argument("--jackknife_min_fraction_passed", help="what fraction of all iterations should"
                                                                "pass to consider the motif robust", type=float)
    parser.add_argument("--job_id", help="job id", type=int)


    parser.set_defaults(
        rna_fastafile='tests/data/test_seqs.fa',
        exp_values_file='tests/data/test_expression_values.csv',
        anno_name_column='ensembl',
        measur_column='measurement',
        seeds_file='tests/data/test_seeds_20.bin',
        temp_folder='tests/temp',
        out='tests/out',
        nbins=5,
        min_occurences=5,
        maxfreq=0.5,
        are_input_seeds_degenerate = False,
        n_permutations=1000,
        max_pvalue=0.01,
        min_zscore=-1,

        step_1_jump=5,
        step_2_min_interval=1,
        step_1_min_fraction=0.8,
        step_2_min_fraction=0.8,
        step_3_min_fraction=0.9,
        min_ratio=5,
        random_noseed=0,
        jackknife_n_samples=10,
        jackknife_fraction_retain=0.66,
        jackknife_n_permutations=1000,
        jackknife_max_pvalue=0.0001,
        jackknife_min_fraction_passed=0.6,

        indices_mode = True,
        index_bit_width = 24,
        job_id = 1,
    )

    args = parser.parse_args()
    return args


def generate_filenames(args):
    # list filenames
    temp_folder = args.temp_folder
    out_folder = args.out
    inputs_folder = os.path.join(temp_folder, "inputs")
    interm_folder = os.path.join(temp_folder, "interm")
    rna_bin_filename = os.path.join(inputs_folder, "rna.bin")
    exp_mask_filename = os.path.join(inputs_folder, "exp_mask.bin")
    profiles_filename = os.path.join(interm_folder, "profiles.bin")
    MI_filename = os.path.join(interm_folder, "MI_values.bin")
    indices_filename = os.path.join(interm_folder, "indices.bin")
    passed_seeds_filename = os.path.join(interm_folder, "passed_seeds_%d.bin" % args.job_id)
    passed_profiles_filename = os.path.join(interm_folder, "passed_profiles_%d.bin" % args.job_id)
    combined_passed_seeds_filename = os.path.join(interm_folder, "combined_passed_seeds.bin")
    combined_passed_profiles_filename = os.path.join(interm_folder, "combined_passed_profiles.bin")
    unique_seeds_filename = os.path.join(interm_folder, "unique_seeds.bin")
    unique_profiles_filename = os.path.join(interm_folder, "unique_profiles.bin")
    families_classification_filename = os.path.join(interm_folder, "families_classification.bin")
    optimized_passed_seeds_filename = os.path.join(interm_folder, "optimized_seeds_%d.bin" % args.job_id)
    optimized_passed_profiles_filename = os.path.join(interm_folder, "optimized_profiles_%d.bin" % args.job_id)
    combined_MI_pv_zscores_filename = os.path.join(interm_folder, "optimized_MI_pv_zscores_%d.bin" % args.job_id)
    combined_robustness_filename = os.path.join(interm_folder, "robustness_array_%d.bin" % args.job_id)
    out_info_filename = os.path.join(out_folder, "pyteiser_info.bin")
    out_matches_filename = os.path.join(out_folder, "pyteiser_matches.bin")

    # make directories
    IO.create_folder(temp_folder)
    IO.create_folder(inputs_folder)
    IO.create_folder(interm_folder)
    IO.create_folder(out_folder)
    with open(indices_filename, 'w') as wf:
        wf.write(str(args.job_id))

    filenames_dict = {
        "rna_bin" : rna_bin_filename,
        "exp_mask" : exp_mask_filename,
        "profiles" : profiles_filename,
        "MI" : MI_filename,
        "indices" : indices_filename,
        "passed_seeds" : passed_seeds_filename,
        "passed_profiles" : passed_profiles_filename,
        "combined_passed_seeds": combined_passed_seeds_filename,
        "combined_passed_profiles": combined_passed_profiles_filename,
        "unique_seeds" : unique_seeds_filename,
        "unique_profiles" : unique_profiles_filename,
        "families_classification" : families_classification_filename,
        "optimized_passed_seeds": optimized_passed_seeds_filename,
        "optimized_passed_profiles": optimized_passed_profiles_filename,
        "combined_MI_pv_zscores" : combined_MI_pv_zscores_filename,
        "combined_robustness" : combined_robustness_filename,
        "out_info" : out_info_filename,
        "out_matches" : out_matches_filename,
        "inputs_folder" : inputs_folder,
        "interm_folder" : interm_folder
    }

    return filenames_dict

def clean_temp_files(filenames_dict):
    for fn in filenames_dict:
        if fn not in ["out_info", "out_matches", "inputs_folder", "interm_folder"]:
            os.remove(filenames_dict[fn])
        elif fn in ["inputs_folder", "interm_folder"]:
            os.rmdir(filenames_dict[fn])

def main():
    args = handler()
    filenames_dict = generate_filenames(args)

    binarize_sequences_args = ["--rna_fastafile", args.rna_fastafile, "--rna_bin_file", filenames_dict['rna_bin']]
    preprocess_custom_expression_profile_args = ["--rna_bin_file", filenames_dict['rna_bin'], "--exp_values_file",
        args.exp_values_file, "--exp_mask_file", filenames_dict['exp_mask'], "--anno_name_column",
        args.anno_name_column, "--measur_column", args.measur_column
    ]
    combine_passed_seeds_args = ["--passed_seed_folder", filenames_dict['interm_folder'],
        "--passed_profiles_folder", filenames_dict['interm_folder'], "--passed_seed_filename_template",
        'passed_seeds', "--passed_profiles_filename_template", "passed_profiles",
        "--combined_seeds_filename", filenames_dict['combined_passed_seeds'], "--combined_profiles_filename",
        filenames_dict['combined_passed_profiles'], "--input_indices_list_file", filenames_dict['indices'],
        "--indices_mode", str(args.indices_mode), "--index_bit_width", str(args.index_bit_width)
    ]
    filter_passed_seeds_with_CMI_args = ["--combined_seeds_filename", filenames_dict['combined_passed_seeds'],
        "--combined_profiles_filename", filenames_dict['combined_passed_profiles'], "--unique_seeds_filename",
        filenames_dict['unique_seeds'], "--unique_profiles_filename", filenames_dict['unique_profiles'],
        "--families_classification_filename", filenames_dict['families_classification'],"--exp_mask_file",
        filenames_dict['exp_mask'], "--nbins", str(args.nbins), "--min_ratio", str(args.min_ratio),
        "--indices_mode", str(args.indices_mode), "--index_bit_width", str(args.index_bit_width),
        "--do_print", str(True)
    ]
    report_csvs_args = ["--rna_fastafile", args.rna_fastafile, "--rna_bin_file", filenames_dict['rna_bin'],
        "--exp_mask_file", filenames_dict['exp_mask'], "--combined_seeds_filename", filenames_dict['optimized_passed_seeds'],
        "--combined_profiles_filename", filenames_dict['optimized_passed_profiles'], "--combined_MI_pv_zscores_filename",
        filenames_dict['combined_MI_pv_zscores'], "--combined_robustness_filename",
        filenames_dict['combined_robustness'], "--out_info_table", filenames_dict['out_info'],
        "--out_matches_table", filenames_dict['out_matches'], "--indices_mode", str(args.indices_mode)
    ]

    binarize_sequences.main(binarize_sequences_args)
    preprocess_custom_expression_profile.main(preprocess_custom_expression_profile_args)
    calculate_seed_profiles.non_sge_dependent_main(
                           args.seeds_file, filenames_dict['profiles'], filenames_dict['rna_bin'],
                           args.are_input_seeds_degenerate, args.indices_mode, args.index_bit_width)
    calculate_MI_profiles.non_sge_dependent_main(
                            filenames_dict['profiles'], filenames_dict['MI'],
                            filenames_dict['exp_mask'],
                            args.indices_mode,
                            args.nbins,
                            args.min_occurences,
                            do_print = True)
    choose_significant_seeds_v3.non_sge_dependent_main(args.seeds_file,
                           filenames_dict['profiles'],
                           filenames_dict['MI'],
                           filenames_dict['passed_seeds'],
                           filenames_dict['passed_profiles'],
                           filenames_dict['exp_mask'],
                           args.max_pvalue, args.min_zscore, args.n_permutations,
                           args.step_1_min_fraction, args.step_1_jump, args.step_2_min_interval,
                           args.step_3_min_fraction,
                           args.indices_mode,
                           args.index_bit_width,
                           do_print=True)
    combine_passed_seeds.main(combine_passed_seeds_args)
    filter_passed_seeds_with_CMI.main(filter_passed_seeds_with_CMI_args)
    optimize_seeds_single_chunk.non_sge_dependent_main(
        str(args.job_id),
        filenames_dict['rna_bin'],
        filenames_dict['exp_mask'],
        filenames_dict['unique_seeds'],
        args.nbins,
        filenames_dict['unique_profiles'],
        args.indices_mode,
        args.index_bit_width,
        0,
        "optimized_seeds",
        "optimized_profiles",
        "optimized_MI_pv_zscores",
        "robustness_array",
        filenames_dict['interm_folder'],
        filenames_dict['interm_folder'],
        filenames_dict['interm_folder'],
        filenames_dict['interm_folder'],
        args.min_occurences,
        args.maxfreq,
        args.n_permutations,
        args.random_noseed,
        args.jackknife_n_permutations,
        args.jackknife_max_pvalue,
        args.jackknife_n_samples,
        args.jackknife_fraction_retain,
        args.jackknife_min_fraction_passed,
        do_chunk_seeds = False
    )
    report_csvs.main(report_csvs_args)
    clean_temp_files(filenames_dict)


if __name__ == '__main__':
    main()


