def invariant(f, x, f_i, y):
    assert f(x) == y
    assert f_i(y) == x

def degenerate_equal(f, coordinate, locations):
    results = [f(loc) for loc in locations]

    # First condition: first maps correctly
    assert results[0] == coordinate, (
        f"\nFirst location: {locations[0]}"
        f"\nExpected:       {coordinate}"
        f"\nGot:            {results[0]}"
    )

    # Second condition: all map to same coordinate
    assert len(set(results)) == 1, (
        f"\nLocations: {locations}"
        f"\nResults:   {results}"
        f"\nExpected all to map to the same coordinate"
    )