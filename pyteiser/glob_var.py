import numpy as np

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
STRUCT_LIST_CHAR = ['stem', 'loop']

_left_stem = np.uint(5)
_right_stem = np.uint(6)
_char_to_extended_structure = {
                    "(" : _left_stem,
                    "<" : _left_stem,
                    "." : _loop,
                    ")" : _right_stem,
                    ">" : _right_stem
                    }

_extended_structure_to_char = {
    _left_stem : "<",
    _right_stem : ">",
    _loop : "."
}

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

_degenerate_nts_mapping = {_U : set([_U]),
                       _C : set([_C]),
                       _G : set([_G]),
                       _A : set([_A]),
                       _N : set([]),
                       _Y : set([_U, _C]),
                       _R : set([_A, _G]),
                       _K : set([_U, _G]),
                       _M : set([_A, _C]),
                       _S : set([_G, _C]),
                       _W : set([_A, _U]),
                       _B : set([_G, _U, _C]),
                       _D : set([_G, _A, _U]),
                       _H : set([_A, _C, _U]),
                       _V : set([_G, _C, _A])}


# see Thruston's responce here for set initiation with single string argument https://stackoverflow.com/questions/17373161/use-curly-braces-to-initialize-a-set-in-python
_degenerate_nts_char_mapping = {'U' : set('U'),
                       'T' : set('U'),
                       'C' : set('C'),
                       'G' : set('G'),
                       'A' : set('A'),
                       'N' : set(''),
                       'Y' : set('UC'),
                       'R' : set('AG'),
                       'K' : set('UG'),
                       'M' : set('AC'),
                       'S' : set('GC'),
                       'W' : set('AU'),
                       'B' : set('GUC'),
                       'D' : set('GAU'),
                       'H' : set('ACU'),
                       'V' : set('GCA')}

_complementary_nt_sets_dict = {
                        _U : set([_A, _G]),
                        _G : set([_U, _C]),
                        _C : set([_G]),
                        _A : set([_U])}

_complementary_deg_nt_dict = {
                        _U: _R,
                        _C: _G,
                        _G: _Y,
                        _A: _U,
                        _N: _N,
                        _Y: _R,
                        _R: _Y,
                        _K: _N,
                        _M: _K,
                        _S: _B,
                        _W: _D,
                        _B: _N,
                        _D: _N,
                        _H: _D,
                        _V: _B}


_rna_alphabet_list = ['A','C','G','U']
_rna_alphabet_nt_mapping = {_A: 0,
                            _C: 1,
                            _G: 2,
                            _U: 3}
_rna_alphabet_string = ''.join(_rna_alphabet_list)


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
mask_out_seed_value = np.float64(-1)