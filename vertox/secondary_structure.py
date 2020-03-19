from .linearfold import linearfold
import random
import multiprocessing
from .path import PATH_LINEARFOLD_EXECUTABLE

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


def stochastic_helix_disruption(sequence, helix, total_iterations):
    left_indices = [helix[0][0], helix[-1][0]]
    right_indices = [helix[-1][1], helix[0][1]]

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
    p.map(linearfold, scrambled_candidates, [PATH_LINEARFOLD_EXECUTABLE]*len(scrambled_candidates))
    best_recovery = 0
    best_scrambled_recovery_sequence = None
    for i in scrambled_candidates:
        tmp_structure, tmp_energy = linearfold(i)
        tmp_recovery = compute_bp_recovery(find_pairs(structure), find_pairs(tmp_structure), ignore_helix=helix)
        amount_disrupted = compute_bp_disruption(find_pairs(structure), find_pairs(tmp_structure), helix)
        if amount_disrupted >= disruption_lower_bound and tmp_recovery > best_recovery:
            best_scrambled_recovery_sequence = i
            best_scrambled_recovery_structure = tmp_structure
            best_recovery = tmp_recovery
