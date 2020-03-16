def compute_bp_disruption(pairs1, pairs2, helix):
    """
    Counts up the number of bps that are disrupted
    assumes pairs1 is the reference
    returns fraction of pairs disrupted
    """
    total_bps = 2 * len(helix)
    collapsed = [i for sublist in helix for i in sublist]

    disrupted_bps = 0
    for x in collapsed:
        if pairs2[x] is None:
            disrupted_bps += 1
    return disrupted_bps / (float(total_bps))
