import numpy as np
import os
import sys

# to make sure relative imports work when some of the wrappers is being implemented as a script
# see more detailed explanation in the test files

current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)


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

_nt_to_char_mapping = {'U' : _U,
                       'T' : _U,
                       'C' : _C,
                       'G' : _G,
                       'A' : _A,
                       'N' : _N}

_char_to_struct_mapping = {_stem : '<',
                           _loop : '.'}


_paired_probabilities_dict = {_U : -3.0, # np.log2( 0.25 * 0.5), since U can create a base pair with either A or G
                            _C : -4.0, # np.log2( 0.25 * 0.25), since C can create a base pair with G
                            _G : -3.0, # np.log2( 0.25 * 0.5), since G can create a base pair with either U or C
                            _A : -4.0, # np.log2( 0.25 * 0.25), since A can create a base pair with U
                            _N : -1.415, # np.log2(2*( 0.25 * 0.5 ) + 2*( 0.25 * 0.25 )), since N itself doesn't tell us anything,
                                                                # but the fact that it's paired puts additional constraints,
                                                                # therefore, sum of all the above
                            'loop': np.log2( 0.25 )
                            }

MAX_SEQ_NAME_LENGTH = 50