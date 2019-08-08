import numpy as np

_U = np.uint8(1)
_C = np.uint8(2)
_G = np.uint8(3)
_A = np.uint8(4)
_N = np.uint8(5)

_stem = np.uint8(1)
_loop = np.uint8(2)

_increment = np.uint8(1)

_char_to_nt_mapping = {_U : 'U',
                        _C : 'C',
                        _G : 'G',
                        _A : 'A',
                        _N : 'N'}

_char_to_struct_mapping = {_stem : '<',
                           _loop : '.'}

_paired_probabilities_dict = {_U : np.log2( 0.25 * 0.5), # U can create a base pair with either A or G
                            _C : np.log2( 0.25 * 0.25), # C can create a base pair with G
                            _G : np.log2( 0.25 * 0.5), # G can create a base pair with either U or C
                            _A : np.log2( 0.25 * 0.25), # A can create a base pair with U
                            _N : np.log2(1), # N can interact with anything
                            'loop': np.log2( 0.25 )
                            }