def invariant(f, x, f_i, y):
    args = x if isinstance(x, tuple) else (x,)
    assert f(*args) == y
    assert f_i(y) == args[0]


def degenerate_equal(f, coordinate, locations):
    assert f(locations[0]) == coordinate
    assert len(set(map(f, locations))) == 1
