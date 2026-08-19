def invariant(f, x, f_i, y):
    assert f(x) == y
    assert f_i(y) == x


def degenerate_equal(f, coordinate, locations):
    assert f(locations[0]) == coordinate
    assert len(set(obj.coordinate for obj in map(f, locations))) == 1
