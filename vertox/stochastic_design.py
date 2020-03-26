from .linearfold import linearfold
import random
from .score_structure_perturbations import compute_bp_disruption, compute_bp_recovery, compute_edit_distance
from .secondary_structure import find_pairs, find_helices, generate_variant_dict

import math
import pandas as pd
from p_tqdm import p_map
from itertools import permutations

def stochastic_helix_remodeling(sequence, structure, helix, total_iterations, ignore_helix=None):
    left_indices = [helix[0][0], helix[-1][0]]
    right_indices = [helix[-1][1], helix[0][1]]

    left_original = sequence[left_indices[0]:left_indices[1] + 1]
    right_original = sequence[right_indices[0]:right_indices[1] + 1]

    # If we have too many iterations for a small helix we will effectively try all scrambles
    # This will cause a problem when we generate the scrambles unless we adjust the number of iterations
    # We want to count the number of UNIQUE permutations

    left_max_permutations = len(set(permutations(list(left_original))))
    right_max_permutations = len(set(permutations(list(right_original))))

    maximum_possible_scrambles = left_max_permutations*right_max_permutations
    if total_iterations > maximum_possible_scrambles - 1:
        total_iterations = maximum_possible_scrambles - 1

    # First we generate a bunch of scrambled sequences
    iterations = 0
    scrambled = list(sequence)
    left_sequence = scrambled[left_indices[0]:left_indices[1] + 1]
    right_sequence = scrambled[right_indices[0]:right_indices[1] + 1]
    scrambled_candidates = []

    iterations = 0
    n_scrambles = 0
    while iterations < total_iterations:
        tmp_scrambled = list(sequence)
        random.shuffle(left_sequence)
        random.shuffle(right_sequence)

        tmp_scrambled[left_indices[0]:left_indices[1] + 1] = left_sequence
        tmp_scrambled[right_indices[0]:right_indices[1] + 1] = right_sequence
        tmp_scrambled = ''.join(tmp_scrambled)

        # Only save unique scrambles
        if tmp_scrambled not in scrambled_candidates and tmp_scrambled != sequence:
            scrambled_candidates.append(tmp_scrambled)
            iterations += 1

    # Then we fold each one and see what we get
    # THIS IS THE EXPENSIVE/PARALLELIZED PART
    folded_scrambles = p_map(linearfold, scrambled_candidates)
    scrambled_structures = [x[0] for x in folded_scrambles]
    scrambled_energies = [x[1] for x in folded_scrambles]

    scrambled_results = []
    for scrambled_candidate, fold  in zip(scrambled_candidates, folded_scrambles):
        scrambled_results.append(generate_variant_dict(sequence, structure, scrambled_candidate, fold, helix, ignore_helix=ignore_helix))
    scrambled_results = pd.DataFrame(scrambled_results)

    return scrambled_results


def stochastic_helix_disruption(sequence, structure, helix, total_iterations):
    # Generate a bunch of structures and score them based on how well they disrupt
    scrambled_results = stochastic_helix_remodeling(sequence, structure, helix, total_iterations, ignore_helix=helix)
    scrambled_results['Variant Type'] = 'Stochastic Helix Disruption'
    scrambled_disruption_results = scrambled_results.sort_values(by=['Disruption Score', 'Recovery Score', 'Edit Distance'],
                                                                 ascending=False,
                                                                 ignore_index=True)
    return scrambled_disruption_results

def stochastic_helix_recovery(sequence, structure, helix, total_iterations):
    # Generate a bunch of structures and score them based on how well they disrupt
    scrambled_results = stochastic_helix_remodeling(sequence, structure, helix, total_iterations, ignore_helix=None)
    scrambled_results['Variant Type'] = 'Stochastic Helix Recovery'
    scrambled_disruption_results = scrambled_results.sort_values(by=['Recovery Score', 'Edit Distance'],
                                                                 ascending=False,
                                                                 ignore_index=True)
    return scrambled_disruption_results

