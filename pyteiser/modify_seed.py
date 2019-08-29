import os
import sys
import glob_var









# this function creates an array of motifs having all possible combinations of bases at the specified position
# modified_motifs will have 15 members (including all the degenerate nucleotides)
# source_motif itself will also be one of the members of modified_motifs

def modify_base(source_motif, position):
    # check if position fits falls into the right range
    if position < 0 or position > source_motif.length:
        print("The position you entered falls outside of the motif!")
        sys.exit(1)

    modified_motifs =

    for nt in glob_var.NT_LIST:
