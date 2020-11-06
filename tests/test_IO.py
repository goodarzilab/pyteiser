import argparse
import numpy as np
import timeit

import os
import sys

# to make sure relative import works in order to import test data
current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)
os.chdir(package_home_path)

import pyteiser.glob_var as glob_var
import pyteiser.structures as structures
import pyteiser.IO as IO
import pyteiser.matchmaker as matchmaker



def handler():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rna_fastafile", type=str)
    parser.add_argument("--outfile", type=str)


    parser.set_defaults(
        rna_fastafile='tests/data/test_seqs.fa',
        rna_bin_file='tests/data/test_seqs.bin',
    )

    args, unknown = parser.parse_known_args()

    return args


def test_main():
    args = handler()
    with open(args.rna_bin_file, 'rb') as rb:
        bitstring = rb.read()
        seq_objects_dict, seq_objects_order = IO.decompress_named_sequences(bitstring)
        full_string = IO.write_named_seq_to_fasta(seq_objects_dict, seq_objects_order)
    with open(args.rna_fastafile, 'r') as rf:
        full_fasta_string = rf.read()

    full_fasta_string_Us = full_fasta_string.replace('T','U').replace('ENSU','ENST')
    assert(len(full_string) == len(full_fasta_string_Us))

    assert(full_string == full_fasta_string_Us)


if __name__ == "__main__":
    test_main()
