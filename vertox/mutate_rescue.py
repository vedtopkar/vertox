import subprocess
import random
from Bio.pairwise2 import format_alignment
from .secondary_structure import find_pairs, find_helices, disrupt_helix_flip, recover_helix_flip
from .score_structure_perturbations import compute_bp_disruption, compute_bp_recovery

def generate_mutate_rescue_library(sequence, structure, gc_rescue=True, iterations=500, stochastic_results=1):
    # Generate helix variants
    pairs = find_pairs(structure)
    helices = find_helices(pairs)

    for helix in helices:
        generate_helix_variants(sequence, structure, helix, iterations=iterations, gc_rescue=gc_rescue)


def generate_mutate_rescue_helix(sequence, structure, helix, gc_rescue=True, iterations=500, stochastic_results=1):
    # Generate flipped disruption and rescue sequences
    left_flipped = disrupt_helix_flip(sequence, helix, side='left')
    right_flipped = disrupt_helix_flip(sequence, helix, side='right')
    recovery_helix_flipped = recover_helix_flip(sequence, helix)


    disruption_lower_bound = max(compute_bp_disruption(find_pairs(structure), find_pairs(left_flipped), helix),
                                 compute_bp_disruption(find_pairs(structure), find_pairs(right_flipped), helix))


def generate_helix_flips(sequence, structure, helix, gc_rescue=True, scramble_iterations=500):
    print("Generating helix variants for\nsequence: {}\nstructure:{}\nhelix: {}".format(sequence, structure, helix))


    alignments = pairwise2.align.globalxx(sequence, right_flipped)
    print(format_alignment(*alignments[0]))

    results = {'sequence': sequence,
               'structure': structure,
               'helix': helix,
               'left_WT': left_original,
               'right_WT': right_original,
               'left_flipped': left_flipped,
               'left_flipped_structure': linearfold(left_flipped)[0],
               'right_flipped': right_flipped,
               'right_flipped_structure': linearfold(right_flipped)[0],
               'flipped_rescue': rescued,
               'flipped_rescue_structure': linearfold(rescued)[0]}

    disruption_lower_bound = max(compute_bp_disruption(find_pairs(structure), find_pairs(left_flipped), helix),
                                 compute_bp_disruption(find_pairs(structure), find_pairs(right_flipped), helix))

