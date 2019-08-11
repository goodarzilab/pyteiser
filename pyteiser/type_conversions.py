import structures

def w_to_n_motif(curr_w_motif):
    curr_n_motif = structures.n_motif(curr_w_motif.stem_length,
                                      curr_w_motif.loop_length,
                                      curr_w_motif.sequence,
                                      curr_w_motif.structure)
    return curr_n_motif


def w_to_n_sequence(current_w_sequence):
    curr_n_sequence = structures.n_sequence(current_w_sequence.length,
                                            current_w_sequence.nts)
    return curr_n_sequence


def w_to_n_motifs_list(inp_list):
    new_list = [0] * len(inp_list)
    for ind, curr_w_motif in enumerate(inp_list):
        new_list[ind] = w_to_n_motif(curr_w_motif)
    return new_list


def w_to_n_sequences_list(inp_list):
    new_list = [0] * len(inp_list)
    for ind, curr_w_motif in enumerate(inp_list):
        new_list[ind] = w_to_n_sequence(curr_w_motif)
    return new_list