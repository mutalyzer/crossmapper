def invariant(f, x, f_i, y):
    assert f(x) == y
    assert f_i(y) == x


def degenerate_equal(f, coordinate, locations):
    assert f(locations[0]) == coordinate
    assert len(
        set(map(f, locations))) == 1


def raise_error(f, x, msg):
    try:
        f(x)
    except ValueError as error:
        assert error.args[0] == msg
        return
    raise ValueError('no error raised')
