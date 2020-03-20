from vertox.linearfold import linearfold
from multiprocessing import Pool
from vertox.path import PATH_LINEARFOLD_EXECUTABLE
from vertox.secondary_structure import *
import pandas as pd

sequence = 'TTCTAATACGACTCACTATAGGCCAAAGGCGTCGAGTAGACGCCAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCGTCAGCGAGTAGCTGACAAAAAGAAACAACAACAACAAC'
structure ='...........................((((((.....))))))...........((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...)).............((((((.....))))))......................'

pairs = find_pairs(structure)
helices = find_helices(pairs)

helix = helices[-2]

disruption_candidates = stochastic_helix_disruption(sequence, structure, helix, 1000)
disruption_candidates.to_csv('test.csv')
