from vertox.secondary_structure import *
from vertox.stochastic_design import *
from vertox.mutate_rescue import *

sequence = 'TTCTAATACGACTCACTATAGGCCAAAGGCGTCGAGTAGACGCCAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCGTCAGCGAGTAGCTGACAAAAAGAAACAACAACAACAAC'
structure ='...........................((((((.....))))))...........((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...)).............((((((.....))))))......................'

# sequence = 'AAAAAAAA'
# structure = '(((..((..))..)))'
# pairs = find_pairs(structure)
# helices = find_helices(pairs)
#
# print(pairs)
# print(helices)

generate_mutate_rescue_library(sequence, structure, total_iterations=25)
