from .secondary_structure import find_pairs, find_helices, disrupt_helix_flip, recover_helix_flip, generate_variant_dict
from .score_structure_perturbations import compute_bp_disruption, compute_bp_recovery, compute_edit_distance
from .stochastic_design import stochastic_helix_disruption, stochastic_helix_recovery
from .linearfold import linearfold

import pandas as pd

def generate_mutate_rescue_library(sequence, structure, gc_rescue=True, total_iterations=500, stochastic_results=1, min_helix_size=3):
    print('Generating mutate-rescue library for:')
    print(sequence)
    print(structure)
    # Generate helix variants
    pairs = find_pairs(structure)
    helices = find_helices(pairs)

    variants = []
    print(helices)
    for helix in helices:
        if len(helix) >= min_helix_size:
            print('Generating mutate rescue for helix: {}'.format(helix))
            # Generate helix flip variants and add them to the list in the form of a DataFrame
            helix_flip_variants = generate_helix_flip_variants(sequence, structure, helix, gc_rescue=gc_rescue)
            variants.append(pd.DataFrame(helix_flip_variants))

            # Generate the stochastically generated mutate-rescue sequences and add to the list
            stochastic_disruption_candidates = stochastic_helix_disruption(sequence, structure, helix, total_iterations=total_iterations)
            stochastic_disruption_subset = stochastic_disruption_candidates[:stochastic_results]
            stochastic_disruption_subset['Variant Type'] = 'Stochastic Helix Disruption'
            variants.append(stochastic_disruption_subset)

            stochastic_recovery_candidates = stochastic_helix_recovery(sequence, structure, helix, total_iterations=total_iterations)
            stochastic_recovery_subset = stochastic_recovery_candidates[:stochastic_results]
            stochastic_recovery_subset['Variant Type'] = 'Stochastic Helix Rescue'
            variants.append(stochastic_recovery_subset)

    variants = pd.concat(variants)
    variants.to_csv('test_output.csv')

def generate_helix_flip_variants(sequence, structure, helix, gc_rescue=True):
    # First we generate the helix flip disruption and recovery sequences
    left_flipped_sequence = disrupt_helix_flip(sequence, helix, side='left')
    right_flipped_sequence = disrupt_helix_flip(sequence, helix, side='right')
    helix_rescue_sequence = recover_helix_flip(sequence, helix, gc_rescue=gc_rescue)

    # Pop them into a variant dictionary for nice legible output
    left_flipped_dict = generate_variant_dict(sequence, structure, left_flipped_sequence, None, helix)
    right_flipped_dict = generate_variant_dict(sequence, structure, right_flipped_sequence, None, helix)
    helix_flipped_rescue_dict = generate_variant_dict(sequence, structure, helix_rescue_sequence, None, helix, ignore_helix=helix)

    left_flipped_dict['Variant Type'] = 'Left Flip'
    right_flipped_dict['Variant Type'] = 'Right Flip'
    helix_flipped_rescue_dict['Variant Type'] = 'Flipped Rescue'

    # Return a list of the results
    return [left_flipped_dict, right_flipped_dict, helix_flipped_rescue_dict]


