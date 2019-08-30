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
_N = np.uint8(5) # any

_Y = np.uint8(6) # UC
_R = np.uint8(7) # AG
_K = np.uint8(8) # UG
_M = np.uint8(9) # AC
_S = np.uint8(10) # GC
_W = np.uint8(11) # AU

_B = np.uint8(12) # GUC
_D = np.uint8(13) # GAU
_H = np.uint8(14) # ACU
_V = np.uint8(15) # GCA

NT_LIST = np.array([_U, _C, _G, _A, _N,
                    _Y, _R, _K, _M, _S, _W,
                    _B, _D, _H, _V], dtype=np.uint8)

MAX_NT = _V

_stem = np.uint8(1)
_loop = np.uint8(2)

STRUCT_LIST = np.array([_stem, _loop])

_increment = np.uint8(1)

_char_to_nt_mapping = {_U : 'U',
                        _C : 'C',
                        _G : 'G',
                        _A : 'A',
                        _N : 'N',
                        _Y : 'Y',
                        _R : 'R',
                        _K : 'K',
                        _M : 'M',
                        _S : 'S',
                        _W : 'W',
                        _B : 'B',
                        _D : 'D',
                        _H : 'H',
                        _V : 'V'
                       }

_nt_to_char_mapping = {'U' : _U,
                       'T' : _U,
                       'C' : _C,
                       'G' : _G,
                       'A' : _A,
                       'N' : _N,
                       'Y' : _Y,
                       'R' : _R,
                       'K' : _K,
                       'M' : _M,
                       'S' : _S,
                       'W' : _W,
                       'B' : _B,
                       'D' : _D,
                       'H' : _H,
                       'V' : _V}

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