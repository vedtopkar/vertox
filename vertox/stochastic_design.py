from .linearfold import linearfold
import random
import multiprocessing
from .score_structure_perturbations import compute_bp_disruption, compute_bp_recovery
from .secondary_structure import find_pairs, find_helices

import nwalign3 as nw
import pandas as pd
from p_tqdm import p_map

def stochastic_helix_remodeling(sequence, structure, helix, total_iterations):
    left_indices = [helix[0][0], helix[-1][0]]
    right_indices = [helix[-1][1], helix[0][1]]

    left_original = sequence[left_indices[0]:left_indices[1] + 1]
    right_original = sequence[right_indices[0]:right_indices[1] + 1]

    # First we generate a bunch of scrambled sequences
    iterations = 0
    scrambled = list(sequence)
    left_sequence = scrambled[left_indices[0]:left_indices[1] + 1]
    right_sequence = scrambled[right_indices[0]:right_indices[1] + 1]
    scrambled_candidates = []
    iterations = 0
    while iterations < total_iterations:
        tmp_scrambled = list(sequence)
        random.shuffle(left_sequence)
        random.shuffle(right_sequence)
        if left_sequence != left_original and right_sequence != right_original:
            tmp_scrambled[left_indices[0]:left_indices[1] + 1] = left_sequence
            tmp_scrambled[right_indices[0]:right_indices[1] + 1] = right_sequence
            scrambled_candidates.append(''.join(tmp_scrambled))
            iterations += 1

    # Then we fold each one and see what we get
    folded_scrambles = p_map(linearfold, scrambled_candidates)
    scrambled_structures, scrambled_energies = list(zip(*folded_scrambles))



    # Align each variant with the WT sequence and save the alignment scores
    alignment_scores = p_map(nw.score_alignment, [sequence]*len(scrambled_candidates), scrambled_candidates)


    return scrambled_candidates, scrambled_structures, scrambled_energies, alignment_scores


def stochastic_helix_disruption(sequence, structure, helix, total_iterations):
    # Generate a bunch of structures and score them based on how well they disrupt
    scrambled_candidates, scrambled_structures, scrambled_energies, alignment_scores = stochastic_helix_remodeling(sequence,
                                                                                                                   structure,
                                                                                                                   helix,
                                                                                                                   total_iterations)

    # Then we score the disruption of each
    # And also the degree to which the rest of the structure remains intact
    WT_pairs = find_pairs(structure)
    scrambled_disruption_scores = [compute_bp_disruption(WT_pairs, find_pairs(s), helix) for s in scrambled_structures]
    scrambled_recovery_scores = [compute_bp_recovery(WT_pairs, find_pairs(s), ignore_helix=helix) for s in
                                 scrambled_structures]

    # Put results into a nice Pandas DataFrame, sort by disruption then recovery, and return the result
    scrambled_disruption_results = [{'scrambled_disruption_sequence': scramble_seq,
                                     'scrambled_disruption_structure': scramble_struct,
                                     'scrambled_disruption_energy': scramble_energy,
                                     'scrambled_alignment_score': alignment_score,
                                     'scrambled_disruption_score': scramble_disruption,
                                     'scrambled_recovery_score': scramble_recovery}
                                    for scramble_seq, scramble_struct, scramble_energy, alignment_score, scramble_disruption, scramble_recovery
                                    in zip(scrambled_candidates, scrambled_structures, scrambled_energies, alignment_scores,
                                           scrambled_disruption_scores, scrambled_recovery_scores)]

    scrambled_disruption_results = pd.DataFrame(scrambled_disruption_results)
    scrambled_disruption_results.sort_values(by=['scrambled_disruption_score', 'scrambled_recovery_score'],
                                             ascending=False,
                                             ignore_index=True,
                                             inplace=True)
    return scrambled_disruption_results


def stochastic_recovery_design(sequence, structure, helix, total_iterations):
    # Generate a bunch of structures and score them based on how well they disrupt
    scrambled_candidates, scrambled_structures, scrambled_energies, alignment_scores = stochastic_helix_remodeling(sequence,
                                                                                                                   structure,
                                                                                                                   helix,
                                                                                                                   total_iterations)

    # Then we score the disruption of each
    # And also the degree to which the rest of the structure remains intact
    WT_pairs = find_pairs(structure)
    scrambled_recovery_scores = [compute_bp_recovery(WT_pairs, find_pairs(s)) for s in scrambled_structures]

    # Put results into a nice Pandas DataFrame, sort by disruption then recovery, and return the result
    scrambled_recovery_results = [{'scrambled_recovery_sequence': scramble_seq,
                                     'scrambled_recovery_structure': scramble_struct,
                                     'scrambled_recovery_energy': scramble_energy,
                                     'scrambled_recovery_score': scramble_recovery,
                                     'scrambled_alignment_score': alignment_score,
                                   }
                                    for scramble_seq, scramble_struct, scramble_energy, alignment_score, scramble_recovery
                                    in zip(scrambled_candidates, scrambled_structures, scrambled_energies, alignment_scores, scrambled_recovery_scores)]

    scrambled_recovery_results = pd.DataFrame(scrambled_recovery_results)
    scrambled_recovery_results.sort_values(by=['scrambled_recovery_score', 'scrambled_alignment_score'],
                                             ascending=[False, True],
                                             ignore_index=True,
                                             inplace=True)
    return scrambled_recovery_results

