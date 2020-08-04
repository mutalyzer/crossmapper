from mutalyzer_crossmapper import cut_locations, nearest_location
from mutalyzer_crossmapper.location import _loc, _nearest_boundary


_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_adjacent_exons = [(1, 3), (3, 5)]


def test_loc():
    assert _loc(1, 2) == [(1, 2)]
    assert _loc(1, 1) == []
    assert _loc(2, 1) == []


def test_cut_locations():
    assert cut_locations(_adjacent_exons, 2) == ([(1, 2)], [(2, 3), (3, 5)])
    assert cut_locations(_adjacent_exons, 3) == ([(1, 3)], [(3, 5)])
    assert cut_locations(_adjacent_exons, 1) == ([], [(1, 3), (3, 5)])
    assert cut_locations(_adjacent_exons, 5) == ([(1, 3), (3, 5)], [])


def test_nearest_boundary():
    assert _nearest_boundary(10, 20, 14, 0) == 0
    assert _nearest_boundary(10, 20, 14, 1) == 0
    assert _nearest_boundary(10, 20, 15, 0) == 1
    assert _nearest_boundary(10, 20, 15, 1) == 1

    assert _nearest_boundary(10, 19, 14, 0) == 0
    assert _nearest_boundary(10, 19, 14, 1) == 1


def test_nearest_location():
    assert nearest_location(_exons, 6) == 0
    assert nearest_location(_exons, 42) == 3
    assert nearest_location(_exons, 71) == 5

    assert nearest_location(_exons, 0) == 0
    assert nearest_location(_exons, 90) == 5

    assert nearest_location(_exons, 10) == 0
    assert nearest_location(_exons, 11) == 1

    assert nearest_location(_exons, 37) == 2
    assert nearest_location(_exons, 38) == 3


def test_nearest_location_force_right():
    assert nearest_location(_exons, 10, 1) == 0
    assert nearest_location(_exons, 11, 1) == 1
    assert nearest_location(_exons, 37, 1) == 3
    assert nearest_location(_exons, 36, 1) == 2


def test_nearest_location_adjacent():
    assert nearest_location(_adjacent_exons, 3) == 1


def test_nearest_location_middle():
    assert nearest_location([(3, 6), (8, 13)], 6) == 0
    assert nearest_location([(3, 6), (8, 13)], 6, 1) == 0

    assert nearest_location([(3, 6), (8, 13)], 7) == 1
    assert nearest_location([(3, 6), (8, 13)], 7, 1) == 1

    assert nearest_location([(3, 6), (9, 13)], 7) == 0
    assert nearest_location([(3, 6), (9, 13)], 7, 1) == 1
