from .linearfold import linearfold
import random
import multiprocessing
from .score_structure_perturbations import compute_bp_disruption, compute_bp_recovery

import pandas as pd

def find_pairs(structure):
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


def find_helices(pairs):
    helices = []
    current_helix = []
    previous = None
    for i in range(len(pairs)):
        if pairs[i] is not None:
            left = i
            right = pairs[i]
            if left < right:  # We haven't seen this bp yet
                bp = [left, right]
                if previous is None or len(current_helix) == 0:  # our first bp
                    current_helix.append(bp)
                    previous = bp
                elif bp[0] - previous[0] == 1 and bp[1] - previous[1] == -1:  # this is a consecutive bp
                    current_helix.append(bp)
                    previous = bp
                else:  # non-consecutive bp
                    print(current_helix)
                    helices.append(current_helix)
                    current_helix = [bp]
                    previous = bp

        elif len(current_helix) > 1:
            print(current_helix)
            helices.append(current_helix)
            current_helix = []
        else:
            continue
    return helices


def disrupt_helix_flip(sequence, helix, side):
    # Determine the indices to flip
    if side == 'left':
        indices = [helix[0][0], helix[-1][0]]
    elif side == 'right':
        indices = [helix[-1][1], helix[0][1]]

    # Flip the inputted side
    flipped = list(sequence)
    flipped[indices[0]:indices[1] + 1] = flipped[indices[0]:indices[1] + 1 ][::-1]

    return ''.join(flipped)


def recover_helix_flip(sequence, helix, gc_rescue=True):
    rescued = list(sequence)

    # Determine indices to flip
    left_indices = [helix[0][0], helix[-1][0]]
    right_indices = [helix[-1][1], helix[0][1]]

    # Flip both sides
    rescued[left_indices[0]:left_indices[1] + 1] = rescued[left_indices[0]:left_indices[1] + 1 ][::-1]
    rescued[right_indices[0]:right_indices[1] + 1] = rescued[right_indices[0]:right_indices[1] + 1 ][::-1]

    # Convert GUs to CGs
    if gc_rescue:
        for pair in helix:
            left = pair[0]
            right = pair[1]
            if gc_rescue:
                if [rescued[left], rescued[right]] == ['G', 'U']:
                    rescued[right] = 'C'
                elif [rescued[left], rescued[right]] == ['U', 'G']:
                    rescued[left] = 'C'

    return ''.join(rescued)


def stochastic_helix_disruption(sequence, structure, helix, total_iterations):
    """
    :param sequence: a sequence of interest
    :param structure: the structure of the sequence (in dot-bracket notation)
    :param helix: a list of base-pair indices defining the helix to be disrupted
    :param total_iterations: number of scramble candidates to generate and score
    :return: a Pandas DataFrame sorted by the disruption of the target helix and the preservation of the remaining structure
    """
    left_indices = [helix[0][0], helix[-1][0]]
    right_indices = [helix[-1][1], helix[0][1]]

    left_original = sequence[left_indices[0]:left_indices[1]+1]
    right_original = sequence[right_indices[0]:right_indices[1]+1]

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

    # Then we score each one and see what we get
    p = multiprocessing.Pool()
    folded_scrambles = p.map(linearfold, scrambled_candidates)
    scrambled_structures, scrambled_energies = list(zip(*folded_scrambles))

    # Then we score the disruption of each
    # And also the degree to which the rest of the structure remains intact
    WT_pairs = find_pairs(structure)
    scrambled_disruption_scores = [compute_bp_disruption(WT_pairs, find_pairs(s), helix) for s in scrambled_structures]
    scrambled_recovery_scores = [compute_bp_recovery(WT_pairs, find_pairs(s), ignore_helix=helix) for s in scrambled_structures]

    # Put results into a nice Pandas DataFrame, sort by disruption then recovery, and return the result
    scrambled_disruption_results = [{'scrambled_disruption_sequence': scramble_seq,
                                    'scrambled_disruption_structure': scramble_struct,
                                    'scrambled_disruption_energy': scramble_energy,
                                    'scrambled_disruption_score': scramble_disruption,
                                    'scrambled_recovery_score': scramble_recovery}
                                    for scramble_seq, scramble_struct, scramble_energy, scramble_disruption, scramble_recovery
                                    in zip(scrambled_candidates,scrambled_structures, scrambled_energies, scrambled_disruption_scores,scrambled_recovery_scores)]

    scrambled_disruption_results = pd.DataFrame(scrambled_disruption_results)
    return scrambled_disruption_results.sort_values(by=['scrambled_disruption_score', 'scrambled_recovery_score'],
                                                    ascending=False,
                                                    ignore_index=True)

