from .secondary_structure import find_pairs, find_helices, disrupt_helix_flip, recover_helix_flip
from .score_structure_perturbations import compute_bp_disruption, compute_bp_recovery, compute_edit_distance
from .linearfold import linearfold

def generate_variant_dict(sequence, structure, variant_sequence, variant_structure, helix):
    WT_pairs = find_pairs(structure)

    if variant_structure is None:
        variant_structure = linearfold(variant_sequence)

    variant_pairs = find_pairs(variant_sequence)

    return {
        'Sequence': variant_sequence,
        'Structure': variant_structure[0],
        'Folding Energy': variant_structure[1],
        'Disruption Score': compute_bp_disruption(WT_pairs, variant_pairs, helix),
        'Edit Distance': compute_edit_distance(sequence, variant_sequence)
    }


def generate_mutate_rescue_library(sequence, structure, gc_rescue=True, iterations=500, stochastic_results=1):
    # Generate helix variants
    pairs = find_pairs(structure)
    helices = find_helices(pairs)

    variants = []

    for helix in helices:
        # Generate helix flip variants and add them to the list
        helix_flip_variants = generate_helix_flip_variants(sequence, structure, helix, iterations=iterations, gc_rescue=gc_rescue)
        variants.extend(helix_flip_variants)

        # Generate the stochastically generated mutate-rescue sequences and add to the list



def generate_helix_flip_variants(sequence, structure, helix, gc_rescue=True, iterations=500, stochastic_results=1):
    # First we generate the helix flip disruption and recovery sequences
    left_flipped_sequence = disrupt_helix_flip(sequence, helix, side='left')
    right_flipped_sequence = disrupt_helix_flip(sequence, helix, side='right')
    helix_rescue_sequence = recover_helix_flip(sequence, helix)

    # Pop them into a variant dictionary for nice legible output
    left_flipped_dict = generate_variant_dict(sequence, structure, left_flipped_sequence, None, helix)
    right_flipped_dict = generate_variant_dict(sequence, structure, right_flipped_sequence, None, helix)
    helix_flipped_rescue_dict = generate_variant_dict(sequence, structure, helix_rescue_sequence, None, helix)

    pairs = find_pairs(sequence)

    # Add recovery score to the flipped rescue
    helix_flipped_rescue_dict['Recovery Score'] = compute_bp_recovery(pairs, find_pairs(helix_rescue_sequence))

    # Return a list of the results
    return [left_flipped_dict, right_flipped_dict, helix_flipped_rescue_dict]


def generate_helix_stochastic_variants(sequence, structure, helix, gc_rescue=True, scramble_iterations=500):


    disruption_lower_bound = max(compute_bp_disruption(find_pairs(structure), find_pairs(left_flipped), helix),
                                 compute_bp_disruption(find_pairs(structure), find_pairs(right_flipped), helix))

