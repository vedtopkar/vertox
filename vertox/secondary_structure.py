
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


