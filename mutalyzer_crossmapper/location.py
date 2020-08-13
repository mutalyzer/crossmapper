def _nearest_boundary(lb, rb, c, p):
    """Find the boundary nearest to `c`. In case of a draw, the parameter `p`
    decides which one is chosen.

    :arg int lb: Left boundary.
    :arg int rb: Right boundary.
    :arg int c: Coordinate (`lb` <= `c` <= `rb`)).
    :arg int p: Preference in case of a draw: 0: left, 1: right.

    :returns int: Nearest boundary: 0: left, 1: right.
    """
    dl = c - lb + 1
    dr = rb - c

    if dl < dr:
        return 0
    if dl > dr:
        return 1
    return p


def nearest_location(ls, c, p=0):
    """Find the location nearest to `c`. In case of a draw, the parameter `p`
    decides which index is chosen.

    :arg list ls: List of locations.
    :arg int c: Coordinate.
    :arg int p: Preference in case of a draw: 0: left, 1: right.

    :returns int: Nearest location.
    """
    rb = len(ls) - 1
    lb = 0

    while lb <= rb:
        i = (lb + rb) // 2

        if c < ls[i][0]:     # `c` lies before this location.
            rb = i - 1
        elif c >= ls[i][1]:  # `c` lies after this location.
            lb = i + 1
        else:                # `c` lies in this location.
            return i

    if i and c < ls[i][0]:  # `c` lies before this location.
        return i - 1 + _nearest_boundary(ls[i - 1][1], ls[i][0], c, p)
    if i < len(ls) - 1:     # `c` lies after this location.
        return i + _nearest_boundary(ls[i][1], ls[i + 1][0], c, p)

    return i
