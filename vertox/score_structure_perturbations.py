import regex

def compute_bp_disruption(pairs1, pairs2, helix):
    """
    Counts up the number of bps that are disrupted
    assumes pairs1 is the reference
    returns fraction of pairs disrupted
    """
    total_bps = len(helix)
    collapsed = [i for sublist in helix for i in sublist]

    disrupted_bps = 0
    for x in collapsed:
        if pairs2[x] is None:
            disrupted_bps += 1
    return disrupted_bps / (float(2*total_bps))


def compute_bp_recovery(pairs1, pairs2, ignore_helix=None, specific_helix=None, false_positive_penalty=0.5):
    """
    Counts how many base pairs are recovered
    assumes pairs1 is the reference
    """
    assert len(pairs1) == len(pairs2)
    assert type(pairs1) == type(pairs2) == list
    recovered = 0

    if ignore_helix is not None:
        ignore_helix = [i for sublist in ignore_helix for i in sublist]
        for index in ignore_helix:
            pairs1[index] = pairs2[index] = None

    total_bps = float(len([x for x in pairs1 if x is not None]))

    for i in range(len(pairs1)):
        if pairs1[i] == pairs2[i] == None:
            continue
        elif pairs1[i] == pairs2[i]:
            recovered += 1
        elif pairs1[i] is None and pairs2[i] is not None:  # False positive
            recovered -= false_positive_penalty
    return recovered / total_bps

def compute_edit_distance(sequence1, sequence2):
    """
    Given two sequences
    returns the edit distance
    """
    assert len(sequence1) == len(sequence2)

    edit_distance = 0
    for i in range(len(sequence1)):
        if sequence1[i] != sequence2[i]:
            edit_distance += 1

    return edit_distance

