
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
