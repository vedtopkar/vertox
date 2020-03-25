from vertox.secondary_structure import *
from vertox.stochastic_design import *
from vertox.mutate_rescue import *

sequence = 'TTCTAATACGACTCACTATAGGCCAAAGGCGTCGAGTAGACGCCAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCGTCAGCGAGTAGCTGACAAAAAGAAACAACAACAACAAC'
structure ='...........................((((((.....))))))...........((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...)).............((((((.....))))))......................'

pairs = find_pairs(structure)
helices = find_helices(pairs)


helix = helices[-2]

generate_mutate_rescue_library(sequence, structure, iterations=100)
