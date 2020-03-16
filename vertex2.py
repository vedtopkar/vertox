import subprocess
import random
from Bio.pairwise2 import format_alignment 
from Bio import pairwise2

def generate_helix_flips(sequence, structure, helix, gc_rescue=True, scramble_iterations=500):
    print("Generating helix variants for\nsequence: {}\nstructure:{}\nhelix: {}".format(sequence, structure,helix))

    left_indices = [helix[0][0], helix[-1][0]]
    left_original = sequence[left_indices[0]:left_indices[1]+1]
    right_indices = [helix[-1][1], helix[0][1]]
    right_original = sequence[right_indices[0]:right_indices[1]+1]
    
    # Flip the left half of the helix
    left_flipped = list(sequence)
    left_flipped[left_indices[0]:left_indices[1]+1] = left_flipped[left_indices[0]:left_indices[1]+1][::-1]
    left_flipped = ''.join(left_flipped)
    
    alignments = pairwise2.align.globalxx(sequence, left_flipped) 
    print(format_alignment(*alignments[0])) 
    
    # Flip the right half of the helix
    right_flipped = list(sequence)
    right_flipped[right_indices[0]:right_indices[1]+1] = right_flipped[right_indices[0]:right_indices[1]+1][::-1]
    right_flipped = ''.join(right_flipped)
    
    alignments = pairwise2.align.globalxx(sequence, right_flipped) 
    print(format_alignment(*alignments[0]))
    
    # Flip both sides of the helix
    # optionally, GU pairs => CG pairs for extra stability
    rescued = list(sequence)
    rescued[right_indices[0]:right_indices[1]+1] = rescued[right_indices[0]:right_indices[1]+1][::-1]
    rescued[left_indices[0]:left_indices[1]+1] = rescued[left_indices[0]:left_indices[1]+1][::-1]

    for pair in helix:
        left = pair[0]
        right = pair[1]
        if gc_rescue:
            if [rescued[left],rescued[right]] == ['G','U']:
                rescued[left] = 'G'
                rescued[right] = 'C'

            elif [rescued[left],rescued[right]] == ['U','G']:
                rescued[left] = 'C'
                rescued[right] = 'G'
    rescued = ''.join(rescued)
    
    alignments = pairwise2.align.globalxx(sequence, rescued) 
    print(format_alignment(*alignments[0]))
    
    print('left flipped:\n{}\nright flipped:\n{}\nboth flipped\n{}'.format(left_flipped,
                                                                           right_flipped,
                                                                           rescued))
    
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

    print("Disruption fraction lower bound: {}".format(disruption_lower_bound))
    # Generate scrambled disrupted helix
    iterations = 0
    scrambled = list(sequence)
    left_sequence = scrambled[left_indices[0]:left_indices[1]+1]
    right_sequence = scrambled[right_indices[0]:right_indices[1]+1]
    scrambled_candidates = []
    iterations = 0
    while iterations < scramble_iterations:
        tmp_scrambled = list(sequence)
        random.shuffle(left_sequence)
        random.shuffle(right_sequence)
        if left_sequence != left_original and right_sequence != right_original:
            tmp_scrambled[left_indices[0]:left_indices[1]+1] = left_sequence
            tmp_scrambled[right_indices[0]:right_indices[1]+1] = right_sequence
            scrambled_candidates.append(''.join(tmp_scrambled))
            iterations += 1

    best_recovery = 0
    best_scrambled_recovery_sequence = None
    for i in scrambled_candidates:
        tmp_structure, tmp_energy = linearfold(i)
        tmp_recovery = compute_bp_recovery(find_pairs(structure), find_pairs(tmp_structure), ignore_helix=helix)
        amount_disrupted = compute_bp_disruption(find_pairs(structure), find_pairs(tmp_structure), helix)
        if  amount_disrupted >= disruption_lower_bound and tmp_recovery > best_recovery:
            best_scrambled_recovery_sequence = i
            best_scrambled_recovery_structure = tmp_structure
            best_recovery = tmp_recovery

    if best_scrambled_recovery_sequence is not None:
        results['scrambled_disruption_sequence'] = best_scrambled_recovery_sequence
        results['scrambled_disruption_structure'] = best_scrambled_recovery_structure
    
    
    print("Generating scrambled rescue")
    # Generate a scrambled rescued helix (both sides are scrambled independently) to search for an optimal structure recovery
    if len(helix) > 5:
        scrambled = list(sequence)
        left_sequence = scrambled[left_indices[0]:left_indices[1]+1]
        right_sequence = scrambled[right_indices[0]:right_indices[1]+1]
        scrambled_candidates = []

        iterations = 0
        while iterations < scramble_iterations:
            tmp_scrambled = list(sequence)
            random.shuffle(left_sequence)
            random.shuffle(right_sequence)
            if left_sequence != left_original and right_sequence != right_original:
                tmp_scrambled[left_indices[0]:left_indices[1]+1] = left_sequence
                tmp_scrambled[right_indices[0]:right_indices[1]+1] = right_sequence
                scrambled_candidates.append(''.join(tmp_scrambled))
                iterations += 1

        best_recovery = 0
        for i in scrambled_candidates:
            tmp_structure, tmp_energy = linearfold(i)
            tmp_recovery = compute_bp_recovery(find_pairs(structure), find_pairs(tmp_structure))
            if  tmp_recovery > best_recovery:
                best_scrambled_recovery_sequence = i
                best_scrambled_recovery_structure = tmp_structure
                best_recovery = tmp_recovery
                print(tmp_recovery)
                print(linearfold(i)[0])
                alignments = pairwise2.align.globalxx(sequence, i) 
                print(format_alignment(*alignments[0])) 

        results['scrambled_rescue_sequence'] = best_scrambled_recovery_sequence
        results['scrambled_rescue_structure'] = best_scrambled_recovery_structure

    
    return results
