from vertox.linearfold import linearfold
from multiprocessing import Pool
from vertox.path import PATH_LINEARFOLD_EXECUTABLE
from vertox.secondary_structure import *

sequence = 'TTCTAATACGACTCACTATAGGCCAAAGGCGTCGAGTAGACGCCAACAACGGAATTGCGGGAAAGGGGTCAACAGCCGTTCAGTACCAAGTCTCAGGGGAAACTTTGAGATGGCCTTGCAAAGGGTATGGTAATAAGCTGACGGACATGGTCCTAACCACGCAGCCAAGTCCTAAGTCAACAGATCTTCTGTTGATATGGATGCAGTTCAAAACCAAACCGTCAGCGAGTAGCTGACAAAAAGAAACAACAACAACAAC'
structure ='...........................((((((.....))))))...........((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...)).............((((((.....))))))......................'

pairs = find_pairs(structure)
helices = find_helices(pairs)

p = Pool(24)
candidates = [sequence]*1000
p.map(linearfold, candidates)
