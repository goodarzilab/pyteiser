import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

import glob_var
import structures


# this function creates an array of motifs having all possible combinations of bases at the specified position
# modified_motifs will have 15 members (including all the degenerate nucleotides)
# source_motif itself will also be one of the members of modified_motifs
def modify_base(source_motif, position):
    # check if position fits falls into the right range
    if position < 0 or position > source_motif.length:
        print("The position you entered falls outside of the motif!")
        sys.exit(1)

    modified_motifs = [0] * len(glob_var.NT_LIST)

    for nt in glob_var.NT_LIST:
        curr_change = source_motif.copy()
        curr_change.sequence[position] = nt

    return modified_motifs


# I take advantage of the fact that the only role the stem_length and loop_length variables
# are playing in w_motif class is to initiate an instance with a pre-defined structure.
# Once an instance has been initiated, you can change the structure in any way you would like,
# as long as you make sure that the matchmaker works with such structure properly
def create_template_elongated_motif(source_motif):
    new_stem_length = source_motif.stem_length + 1
    template = structures.w_motif(new_stem_length, source_motif.loop_length)
    template.sequence[0] = glob_var._N
    template.sequence[1:] = source_motif.sequence
    template.structure[0] = glob_var._stem
    template.structure[1:] = source_motif.structure
    return template



# this function creates an array of motifs each having one additional phrase compared to source_motif
# modified_motifs will have 46 members - why??
# modified_motifs will have 1 + 2 * 15 members: the source motif + 2 structures (loop and stem) by
# 15 degenerate nucleotides
# source_motif itself will also be one of the members of modified_motifs
def elongate_motif(source_motif):


    modified_motifs = [0] * 45

    motif_index = 0
    modified_motifs[motif_index] = source_motif.copy()

    for struct_mode in range(1, 3+1):
        for nt in glob_var.NT_LIST:
            motif_index += 1

            structures.w_motif


            # copy source motif
            # increase length by 1
            # assign sequence and structure to the leftmost element
            # maybe create a new class for extended motif?



            #modified_motifs[motif_index] =

