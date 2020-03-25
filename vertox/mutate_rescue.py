from .secondary_structure import find_pairs, find_helices, disrupt_helix_flip, recover_helix_flip
from .score_structure_perturbations import compute_bp_disruption, compute_bp_recovery, compute_edit_distance
from .linearfold import linearfold

class SequenceVariant:
    def __init__(self, WT_sequence, WT_structure):


def generate_variant_dict(sequence, structure, variant_sequence,
                          variant_structure=None, variant_type=None):
    # Type can be mutate or rescue
    assert variant_type is not None

    WT_pairs = find_pairs(sequence)

    if variant_structure is None:
        variant_structure = linearfold(right)

    variant_pairs = find_pairs(variant_sequence)

    return {
        'Sequence': variant_sequence,
        'Structure': variant_structure[0],
        'Folding Energy': variant_folded[1],
        'Disruption Score': compute_bp_disruption(WT_pairs, variant_pairs),
        'Recovery Score': None,
        'Edit Distance': compute_edit_distance(sequence, left_flipped)
    }

def generate_mutate_rescue_library(sequence, structure, gc_rescue=True, iterations=500, stochastic_results=1):
    # Generate helix variants
    pairs = find_pairs(structure)
    helices = find_helices(pairs)

    variants = []

    for helix in helices:
        left_flipped = disrupt_helix_flip(sequence, helix, side='left')
        left_flipped_pairs = find_pairs(left_flipped)
        left_flipped_folded = linearfold(left_flipped)

        right_flipped = disrupt_helix_flip(sequence, helix, side='right')
        right_flipped_pairs = find_pairs(right_flipped)
        right_flipped_folded = linearfold(right_flipped)

        variants.append({
            'Variant Type': 'Left Flipped',
            'Sequence': left_flipped,
            'Structure': left_flipped_folded[0],
            'Folding Energy': left_flipped_folded[1],
            'Disruption Score': compute_bp_disruption(pairs, left_flipped_pairs),
            'Recovery Score': None,
            'Edit Distance': compute_edit_distance(sequence, left_flipped)
        })

        variants.append({
            'Variant Type': 'Right Flipped',
            'Sequence': right_flipped,
            'Structure': right_flipped_folded[0],
            'Folding Energy': right_flipped_folded[1],
            'Disruption Score': compute_bp_disruption(pairs, right_flipped_pairs),
            'Recovery Score': None,
            'Edit Distance': compute_edit_distance(sequence, right_flipped)
        })


        generate_helix_flip_variants(sequence, structure, helix, iterations=iterations, gc_rescue=gc_rescue)



def generate_helix_flip_variants(sequence, structure, helix, gc_rescue=True, iterations=500, stochastic_results=1):
    # Generate flipped disruption and rescue sequences
    left_flipped = disrupt_helix_flip(sequence, helix, side='left')
    right_flipped = disrupt_helix_flip(sequence, helix, side='right')
    recovery_helix_flipped = recover_helix_flip(sequence, helix)

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

    return results


def generate_helix_flips(sequence, structure, helix, gc_rescue=True, scramble_iterations=500):
    print("Generating helix variants for\nsequence: {}\nstructure:{}\nhelix: {}".format(sequence, structure, helix))

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

