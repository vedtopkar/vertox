from .secondary_structure import find_pairs, find_helices
import numpy as np

class SeqStruct:
    def __init__(self, sequence, structure=None, reactivity=None):
        # First, we populate the sequence with the inputted attributes
        self.sequence = sequence
        self.length = len(self.sequence)

        self.structure = structure if structure is not None else '.'*self.length
        self.reactivity = reactivity if reactivity is not None else np.array([np.nan]*self.length)

        # At this point we want to auto-populate te pairs and helices attributes
        self.pairs = find_pairs(self.structure)
        self.helices = find_helices(self.pairs)

        print('Initiated SeqStruct object with sequence: {}\nstructure: {}\nreactivity: {}\n'.format(
            self.sequence,
            self.structure,
            self.reactivity))
        print(self.pairs)
        print(self.helices)

