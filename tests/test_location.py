from mutalyzer_crossmapper import cut_locations, nearest_location
from mutalyzer_crossmapper.location import _loc, _nearest_boundary


def test_loc():
    """Proper locations should remain unaltered, others are removed."""
    assert _loc(1, 2) == [(1, 2)]
    assert _loc(1, 1) == []
    assert _loc(2, 1) == []


def test_cut_locations():
    """A cut should result in two parts containing proper locations."""
    locations = [(1, 3), (3, 5)]

    assert cut_locations(locations, 2) == ([(1, 2)], [(2, 3), (3, 5)])
    assert cut_locations(locations, 3) == ([(1, 3)], [(3, 5)])
    assert cut_locations(locations, 1) == ([], [(1, 3), (3, 5)])
    assert cut_locations(locations, 5) == ([(1, 3), (3, 5)], [])


def test_nearest_boundary_even():
    """Boundaries are not equally near, preference is irrelevant."""
    assert _nearest_boundary(10, 20, 14, 0) == 0
    assert _nearest_boundary(10, 20, 14, 1) == 0
    assert _nearest_boundary(10, 20, 15, 0) == 1
    assert _nearest_boundary(10, 20, 15, 1) == 1


def test_nearest_boundary_odd():
    """Boundaries are equally near, preference is relevant."""
    assert _nearest_boundary(10, 19, 14, 0) == 0
    assert _nearest_boundary(10, 19, 14, 1) == 1


def test_nearest_location():
    """Index of the nearest location."""
    locations = [(10, 20), (30, 40), (50, 60)]

    assert nearest_location(locations, 8) == 0
    assert nearest_location(locations, 15) == 0
    assert nearest_location(locations, 22) == 0

    assert nearest_location(locations, 28) == 1
    assert nearest_location(locations, 35) == 1
    assert nearest_location(locations, 42) == 1

    assert nearest_location(locations, 48) == 2
    assert nearest_location(locations, 55) == 2
    assert nearest_location(locations, 62) == 2


def test_nearest_location_even():
    """Index of the nearest location, preference is irrelevant."""
    assert nearest_location([(3, 6), (8, 13)], 6, 0) == 0
    assert nearest_location([(3, 6), (8, 13)], 6, 1) == 0
    assert nearest_location([(3, 6), (8, 13)], 7, 0) == 1
    assert nearest_location([(3, 6), (8, 13)], 7, 1) == 1


def test_nearest_location_even():
    """Index of the nearest location, preference is relevant."""
    assert nearest_location([(3, 6), (9, 13)], 7) == 0
    assert nearest_location([(3, 6), (9, 13)], 7, 1) == 1


def test_nearest_location_adjacent():
    """Adjacent locations have no overlap."""
    locations = [(1, 3), (3, 5)]

    assert nearest_location(locations, 2) == 0
    assert nearest_location(locations, 3) == 1
