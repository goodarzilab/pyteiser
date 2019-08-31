import os
import sys
import numpy as np
import numba
import argparse

# when you run a python file as a script and you do it from some other folder, relative import simply doesn't work
# a simple explanations can be found here: https://stackoverflow.com/questions/30669474/beyond-top-level-package-error-in-relative-import
# and here https://stackoverflow.com/questions/14132789/relative-imports-for-the-billionth-time/14132912#14132912
# basically, Python discards the knowledge about where the currently running script was imported from and therefore
# relative imports in this script won't work. Again: scripts can't import relative!
# to solve this problem, most people use sys.path hacks. For example, here: https://stackoverflow.com/questions/6323860/sibling-package-imports
# the common sys.path.insert(0, os.path.abspath('..')) trick is not what I want to use, since it adds the folder
# that is one level up from the current working directory, which might be different from where the script you're running is
# what I do instead is I add (to sys.path) the folder one level up from the script I am currently running so it can definitely
# do a relative import
# to find the path of python script currently being executed, see https://stackoverflow.com/questions/595305/how-do-i-get-the-path-of-the-python-script-i-am-running-in
# to go up one folder, see https://stackoverflow.com/questions/9856683/using-pythons-os-path-how-do-i-go-up-one-directory

current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)

import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.matchmaker as matchmaker
import pyteiser.type_conversions as type_conversions
import pyteiser.wrappers.calculate_seed_profiles as calculate_seed_profiles




def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seeds_bin_file", type=str)
    parser.add_argument("--profiles_bin_file", type=str)
    parser.add_argument("--rna_bin_file", type=str)


    parser.set_defaults(
        seeds_bin_file = '/Users/student/Documents/hani/temp/seeds_temp/test_seeds/test_seeds.bin',
        profiles_bin_file='/Users/student/Documents/hani/temp/profiles_temp/test_profiles/test_profiles.bin',
        rna_bin_file='/Users/student/Documents/hani/iTEISER/step_2_preprocessing/reference_files/reference_transcriptomes/binarized/Gencode_v28_GTEx_expressed_transcripts_from_coding_genes_3_utrs_fasta.bin'
    )

    args = parser.parse_args()

    return args




def test_matchmaker():
    # test matchmaking algorithms
    # 3 strings listed here contain instances of 3 matches that are also listed here

    test_motif_1 = structures.w_motif(4,6)
    test_motif_2 = structures.w_motif(4,6)
    test_motif_3 = structures.w_motif(4,6)
    test_motif_1.from_string("GNCANCNNUU")
    test_motif_2.from_string("AAUNNGNGNU")
    test_motif_3.from_string("NNACGNNCUU")
    test_motifs_list_w = [test_motif_1, test_motif_2, test_motif_3]
    test_motifs_list = type_conversions.w_to_n_motifs_list(test_motifs_list_w)

    test_string_1 = 'UUUUUUUGACAACAAUUTGTCUUUUU' # instance motif_1 at 7
    test_string_2 = "GGCAUCAGUUUUUUAAUGUGUGAUCAUUGGGUUCCCCCUUUUU" # instance motif_2 at 14
    test_string_3 = "AAUUAAAACCCCCCCAAACGCCCUUGUUUCCCACCACGGGCUUGUGGAAAAUUUUUU" # instances motif_3 at 15 and 33

    test_sequence_1 = structures.w_sequence(len(test_string_1))
    test_sequence_2 = structures.w_sequence(len(test_string_2))
    test_sequence_3 = structures.w_sequence(len(test_string_3))
    test_sequence_1.from_sequence(test_string_1)
    test_sequence_2.from_sequence(test_string_2)
    test_sequence_3.from_sequence(test_string_3)
    test_sequences_list_w = [test_sequence_1, test_sequence_2, test_sequence_3]
    test_sequences_list = type_conversions.w_to_n_sequences_list(test_sequences_list_w)

    boolean_matchmaker_desired = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=bool)
    boolean_matchmaker_res = np.zeros(shape=(3, 3), dtype=bool)
    indices_matchmaker_desired = [[7],[],[],[],[14],[],[],[],[15,33]]
    indices_matchmaker_res = []

    for i, mt in enumerate(test_motifs_list):
        for k, sq in enumerate(test_sequences_list):
            is_match = matchmaker.is_there_motif_instance(mt, sq)
            matching_indices = matchmaker.find_all_motif_instances(mt, sq)
            boolean_matchmaker_res[i,k] = is_match
            indices_matchmaker_res.append(matching_indices)


    assert(np.array_equal(boolean_matchmaker_res, boolean_matchmaker_desired))
    assert(indices_matchmaker_res == indices_matchmaker_desired)

def test_current_pair(stem = 4, loop = 7,
                      motif_str = "NGCAUNGNANN",
                      seq_str = "UGCAUUGUAUGUGUG"):
    test_motif = structures.w_motif(stem, loop)
    test_motif.from_string(motif_str)
    n_test_motif = type_conversions.w_to_n_motif(test_motif)

    test_sequence = structures.w_sequence(len(seq_str))
    test_sequence.from_sequence(seq_str)
    n_test_sequence = type_conversions.w_to_n_sequence(test_sequence)
    is_match = matchmaker.is_there_motif_instance(n_test_motif, n_test_sequence)

    if is_match:
        print("Sequence %s matches the motif %s" % (seq_str, motif_str))
        motif_instances = matchmaker.find_all_motif_instances(n_test_motif, n_test_sequence)
        print("Motif instances are: ", ", ".join([str(x) for x in motif_instances]))

    else:
        print("Sequence %s DOES NOT matches the motif %s" % (seq_str, motif_str))




def define_constants(args):
    seqs_shape = (4,7)

    seqs_to_test = ['NGCAUNGNANN', 'NACAUNGNANN', 'UNCAUNGNANN', 'CNCAUNGNANN', 'GNCAUNGNANN', 'ANCAUNGNANN', 'NNCAUNGNANN',
     'NUGAUNGNANN', 'NCGAUNGNANN', 'NGGAUNGNANN', 'NAGAUNGNANN', 'UNGAUNGNANN', 'CNGAUNGNANN', 'GNGAUNGNANN',
     'ANGAUNGNANN', 'NNGAUNGNANN', 'NUAAUNGNANN', 'NCAAUNGNANN', 'NGAAUNGNANN', 'NAAAUNGNANN', 'UNAAUNGNANN',
     'CNAAUNGNANN', 'GNAAUNGNANN', 'ANAAUNGNANN', 'NNAAUNGNANN', 'UUNAUNGNANN', 'CUNAUNGNANN', 'GUNAUNGNANN',
     'AUNAUNGNANN', 'NUNAUNGNANN', 'UCNAUNGNANN', 'CCNAUNGNANN', 'GCNAUNGNANN', 'ACNAUNGNANN', 'NCNAUNGNANN',
     'UGNAUNGNANN', 'CGNAUNGNANN', 'GGNAUNGNANN', 'AGNAUNGNANN', 'NGNAUNGNANN', 'UANAUNGNANN', 'CANAUNGNANN',
     'GANAUNGNANN', 'AANAUNGNANN', 'NANAUNGNANN', 'UNNAUNGNANN', 'CNNAUNGNANN', 'GNNAUNGNANN', 'ANNAUNGNANN',
     'NNNAUNGNANN', 'UUUNUNGNANN', 'CUUNUNGNANN', 'GUUNUNGNANN', 'AUUNUNGNANN', 'NUUNUNGNANN', 'UCUNUNGNANN',
     'CCUNUNGNANN', 'GCUNUNGNANN', 'ACUNUNGNANN', 'NCUNUNGNANN', 'UGUNUNGNANN', 'CGUNUNGNANN', 'GGUNUNGNANN',
     'AGUNUNGNANN', 'NGUNUNGNANN', 'UAUNUNGNANN', 'CAUNUNGNANN', 'GAUNUNGNANN', 'AAUNUNGNANN', 'NAUNUNGNANN',
     'UNUNUNGNANN', 'CNUNUNGNANN', 'GNUNUNGNANN', 'ANUNUNGNANN', 'UUCNUNGNANN', 'CUCNUNGNANN', 'GUCNUNGNANN',
     'AUCNUNGNANN', 'NUCNUNGNANN', 'UCCNUNGNANN', 'CCCNUNGNANN', 'GCCNUNGNANN', 'ACCNUNGNANN', 'NCCNUNGNANN',
     'UGCNUNGNANN', 'CGCNUNGNANN', 'GGCNUNGNANN', 'AGCNUNGNANN', 'NGCNUNGNANN', 'UACNUNGNANN', 'CACNUNGNANN',
     'GACNUNGNANN', 'AACNUNGNANN', 'NACNUNGNANN', 'UNCNUNGNANN', 'CNCNUNGNANN', 'GNCNUNGNANN', 'ANCNUNGNANN',
     'NNCNUNGNANN', 'UUGNUNGNANN', 'CUGNUNGNANN']

    bin_file_to_test = args.rna_bin_file

    desired_numbers = [207, 180, 299, 94, 206, 158, 748, 388, 13, 363, 348, 351, 132, 358, 274, 1097, 316, 93, 360, 556, 414, 82, 374, 448,
     1276, 743, 158, 355, 317, 1515, 96, 51, 110, 62, 317, 488, 20, 399, 393, 1260, 443, 202, 368, 481, 1448, 1697, 422,
     1208, 1220, 4146, 2235, 672, 769, 791, 4139, 577, 295, 480, 228, 1552, 1049, 85, 608, 510, 2168, 710, 323, 403,
     547, 1890, 4206, 1353, 2155, 1984, 529, 204, 267, 128, 1097, 242, 132, 191, 80, 643, 272, 24, 363, 181, 828, 151,
     102, 107, 106, 461, 1169, 458, 910, 492, 2875, 850, 473]

    return seqs_shape, seqs_to_test, bin_file_to_test, desired_numbers



def prepare_known_seeds(args):
    seqs_shape, seqs_to_test, bin_file_to_test, desired_numbers = define_constants(args)

    w_motifs_list = [0] * len(seqs_to_test)
    for ind, seq in enumerate(seqs_to_test):
        curr_test_motif = structures.w_motif(seqs_shape[0], seqs_shape[1])
        curr_test_motif.from_string(seq)
        w_motifs_list[ind] = curr_test_motif

    seqs_dict, seqs_order = IO.read_rna_bin_file(bin_file_to_test)
    w_seqs_list = [seqs_dict[name] for name in seqs_order]

    n_motifs_list = type_conversions.w_to_n_motifs_list(w_motifs_list)
    n_seqs_list = type_conversions.w_to_n_sequences_list(w_seqs_list)
    return n_motifs_list, n_seqs_list


def test_calculate_seed_profiles():
    n_motifs_list, n_seqs_list = prepare_known_seeds(args)
    compress_test_seeds(args, n_motifs_list)
    matchmaker.calculate_profiles_list_motifs(n_motifs_list, n_seqs_list, do_print = True)




def compress_test_seeds(args, n_motifs_list):
    motifs_bytestrings_list = [0] * len(n_motifs_list)
    for i, motif in enumerate(n_motifs_list):
        w_motif = type_conversions.n_to_w_motif(motif)
        w_motif.compress()
        motifs_bytestrings_list[i] = w_motif.bytestring
    with open(args.seeds_bin_file, 'wb') as wf:
        string_to_write = b''.join(motifs_bytestrings_list)
        wf.write(string_to_write)




def test_profiles_compression_decompression(args, do_shorten_test = True):
    args.seedfile = args.seeds_bin_file
    n_motifs_list, n_seqs_list = calculate_seed_profiles.prepare_lists_for_calculations(args)
    if do_shorten_test:
        n_motifs_list = n_motifs_list[0:3]
    calculated_profiles_array = calculate_seed_profiles.calculate_write_profiles(n_motifs_list, n_seqs_list,
                                            args.profiles_bin_file, do_print=True,
                                            do_return=True)

    with open(args.profiles_bin_file, 'rb') as rf:
        bitstring = rf.read()
    decompressed_profiles_array = IO.decompress_profiles(bitstring)

    assert(np.array_equal(calculated_profiles_array, decompressed_profiles_array))

    # for i in range(decompressed_profiles_array.shape[0]):
    #     print("Number of mathces in the %d-th array" % i, decompressed_profiles_array[i].sum())
    # for i in range(calculated_profiles_array.shape[0]):
    #     print("Number of mathces in the %d-th array" % i, calculated_profiles_array[i].sum())


if __name__ == "__main__":
    args = handler()

    # test individual cases
    test_matchmaker()

    # test that the number of instances identified is correct
    test_calculate_seed_profiles()

    # test that I can
    test_profiles_compression_decompression(args)



