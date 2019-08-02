class s_phrase:

    def __init__(self, dict):
        self.base = 0 # NUCBIT - byte
        # _U | _C | _G | _A; if structure = _pair, base shows the left nucleotide in the paired bases
        self.structure = 0 # byte: _pair ^ _leftBulge ^ _rightBulge


class s_sequence:

    def __init__(self, dict):
        self.name = ""
        self.bases = "" #list of NUCBITs
        self.length = 0


class s_motif:

    def __init__(self, dict):
        self.phrases = [] # list of s_phrase objects
        self.num_phrases = 0
        self.linear_length = 0



