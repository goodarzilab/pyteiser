import numpy as np

_U = np.uint8(1)
_C = np.uint8(2)
_G = np.uint8(3)
_A = np.uint8(4)
_N = np.uint8(5)

_pair = np.uint8(1)
_leftBulge = np.uint8(2)
_rightBulge = np.uint8(3)

_increment = np.uint8(1)

_char_to_nt_mapping = {_U : 'U',
                        _C : 'C',
                        _G : 'G',
                        _A : 'A',
                        _N : 'N'}

_char_to_struct_mapping = {_pair : '',
                           _leftBulge : '<',
                           _rightBulge : '>'}
