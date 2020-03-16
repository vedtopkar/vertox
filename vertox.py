import numpy as np

class SeqStruct:
    def __init__(self, sequence, structure=None, reactivity=None):
        self.sequence = sequence
        self.length = len(self.sequence)

        if structure is None:
            self.structure = '.'*self.length
        else:
            self.structure = structure

        if reactivity is None:
            self.reactivity = np.array([np.nan]*self.length)
        else:
            self.structure = structure

        self.pairs = self.find_pairs(self.structure)
        self.helices = self.find_helices(self.pairs)
        self.flipped_helices = self.generate_helix_flips(self.sequence, self.helices)

        print('Initiated SeqStruct object with sequence: {}\nstructure: {}\nreactivity: {}\n'.format(
                                                                                                self.sequence,
                                                                                                self.structure,
                                                                                                self.reactivity))
        print(self.pairs)
        print(self.helices)
        print(self.flipped_helices)

    def find_pairs(self, structure):
        stack = []
        pairs = [None]*len(structure)
        for i, c in enumerate(structure):
            if c == '(':
                stack.append(i)
            elif c == ')':
                left = stack.pop()
                right = i
                pairs[right] = left
                pairs[left] = right
            elif c == '.':
                pass
        return pairs

    def find_helices(self, pairs):
        helices = []
        current_helix = []
        for i in range(len(pairs)):
            if pairs[i] is not None:
                left = pairs[i]
                right = pairs[left]
                if left > right:
                    current_helix.append([pairs[i], pairs[pairs[i]]])
            elif len(current_helix) > 1:
                helices.append(current_helix)
                current_helix = []
            else:
                continue
        return helices

    def generate_helix_flips(self, seq, helices):
        flipped_helices = []
        for helix in helices:
            flipped = list(seq)
            for pair in helix:
                i,j = pair
                flipped[i], flipped[j] = flipped[j], flipped[i]
            flipped_helices.append(''.join(flipped))
        return flipped_helices


test = SeqStruct('GAGCCUCCCUGCUCAGCCUUCCCGAAUCCUGCCCUCGGCUUCUUAAUAUAACUGCCUUAAACGUUUAAUUCUACUUGCACCAAAUAGCUAGUUAGAGCAGACCCUCUCUUAAUCCCGUGGGGCUGUGAACGCGGCGGGGCCAGGCCCACGGCACCCUGACUGGCUAAAACUGUUUGUCCCUUUUUAUUUGAAGAUUGAGUUUCCUCGGGGUCUUCUCUGCCCCGACUUGCUCCCCGUGUACCUUGGUCGACUCCGGAGGUUCAGGUGCACGGACACCCUUUCAAGUUCACCCCUACUCCAUCCUCAGACUUUCUUUUCACGGCGAGGCGCACCCCUCCAGCUUCCGUGGGCACUGCGGAUAGACAGGCACACCGCCAAGGAGCCAGAGAGCAUGGCGCAGGGGACUGUGUGGUCCAGGCUUCCUUUGUUUUCUUUCCCCUAAAGAGCUUUGUUUUUCCUAACAGGAUCAGACAGUCUUGGAGUGGCUUACACAACGGGGGCUUGUGGUAUGUGAGCACAGGCUGGGCAGCUGUGAGAGUCCAGAGUGGGGUGGCCCUGGGGACGCUUCCAGGCCAGCGGUUCCCUGCACCCCACCAGCUGAUUUCGAGCGUGGCAGAGGGAAGGAAAGGGGCGAGCGGGCUGGGCAAUGGACCCGACAGGAAACGGGGACUUAGGGGAACACGCUGGAGAUGCCAUGUGUGGCUGCCGAAGGUCACCAUCUCUCCUCAGUGGCUCCCCAGAGCAGGUGCUUUUAAGAACCCUGUUUCCUCUCAGAGCCCAGGGAGAGUCCAAGGACAUGGCGCAUCAGGAAGUGGGACUGCAGGAGUUCUCUGGUGGCCUCGUGCUGUCCCUCUGGCCACUUCUCACUUUAGGGUGGUCAGCGGCAGCUCGCCAUGGCAGUGCCCAUUGGUGCACACUAACCUCAGUGGAAAAGUAACCAUUCCCUGCCUCUUAGAAAGAACUCAUUCUUAGUUUUAGGAGGGUUCCUGUCGCUGAAUCAAGUCGCUGCCCUGGAUGCAGGGCUGGCCUGGGCGACCCUCCAGGGAUGAGGAGCUCAGAAUUCCAGUCUUCUAAUGUCCACGGACACCUCCCCAUCCCUCUAACGUACUGACUAUGUCUUUUGAUUUAGCAUGUCUUCUAUAGACCUUCCAAAGAGACCCACACUGGCACUGUCACCCCCUAGGAGGGAAGGUGAUGGUUGAUGUAGCCCGACGCGCAUCUUGUUAAUCCGUUCUAAUUCCGAGGAGAGUGUGGGUUUAAGAUAACACCUAUUAAUGCAUUGCCACAAUAAUGUGGGGGUAAGAGAAACGCAGGGACGAAACUUCCAGAAACAAACCCUCCAGAUCGUUCCACAGGAGUGUUCGCCCUCCGGUGUGACUGAACGACCGACCUUGCCCAUGGCUUCAUCCAGACAGCACAGCUGCAGUAUGGCUGGACAGAAGCACCUACUGUUCUUGGAUAUUGAAAUAAAAUAAUAAACUUGCAAUGAUCUUU'
                ,'...................................................(((((..............(((((.((........)).)).))).((((.((((.............))))))))......)))))....................(((.(((.......))).))).((((......))))...((((....(((((((.......)))))))...)))).((((((((((.(((((....)))....))))))))))))................................................(((((((......))))..))).......................(((.....)))........................................................................((((((.((((...))))..)))))).......(((.....))).......(((((((.........))))))).......................................................................................................................................(((((((((((.....(((((.(((...(((((....(((.((((((......))))))))).)))))...)))..))))).((((....)))).......))))))))))).................................((((((.....((((((.....)))))))))))).......((((((...(((((((((.((......))))))))))).))))))...................................(((.((((.........))))))).........................................................(((((((((....)))))).))).................................................................................((.(((.(((....))).))))).((((.....)))).........(((((((((((...(((((((((((((...))))..))))))))).(((((.((.....))))))).......(((.............))).)))))))))))........................((((......))))(((((((........((((....................)))).....((((((((((((((........))))..))))...)))))).....)))))))...(((((.(((((.((((......))....)).)))))..)))))...................................................')

