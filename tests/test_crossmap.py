from crossmapper.crossmapper import (
    Crossmap, Locus, MultiLocus, _cut, _loc, _nearest_location, _offsets)


_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_cds = (32, 43)
_adjacent_exons = [(1, 3), (3, 5)]


def _test_invariant(f, x, f_i, y):
    assert f(x) == y
    assert f_i(y) == x


def test_cut():
    assert _cut(_adjacent_exons, 2) == ([(1, 2)], [(2, 3), (3, 5)])
    assert _cut(_adjacent_exons, 3) == ([(1, 3)], [(3, 5)])
    assert _cut(_adjacent_exons, 1) == ([], [(1, 3), (3, 5)])
    assert _cut(_adjacent_exons, 5) == ([(1, 3), (3, 5)], [])


def test_loc():
    assert _loc(1, 2) == [(1, 2)]
    assert _loc(1, 1) == []
    assert _loc(2, 1) == []


def test_offsets():
    assert _offsets(_exons) == [0, 3, 9, 14, 18, 20]
    assert _offsets(_exons, True) == [0, 2, 4, 8, 13, 19]
    assert _offsets(_adjacent_exons) == [0, 2]
    assert _offsets(_adjacent_exons, True) == [0, 2]


def test_nearest_location():
    assert _nearest_location(_exons, 6) == 0
    assert _nearest_location(_exons, 42) == 3
    assert _nearest_location(_exons, 71) == 5

    assert _nearest_location(_exons, 0) == 0
    assert _nearest_location(_exons, 90) == 5

    assert _nearest_location(_exons, 10) == 0
    assert _nearest_location(_exons, 11) == 1
    assert _nearest_location(_exons, 10, 1) == 0
    assert _nearest_location(_exons, 11, 1) == 1

    assert _nearest_location(_exons, 37) == 2
    assert _nearest_location(_exons, 38) == 3
    assert _nearest_location(_exons, 37, 1) == 3
    assert _nearest_location(_exons, 36, 1) == 2

    assert _nearest_location(_adjacent_exons, 3) == 1


def test_Locus():
    locus = Locus(_exons[2])

    _test_invariant(locus.to_position, 29, locus.to_coordinate, (1, -1))
    _test_invariant(locus.to_position, 30, locus.to_coordinate, (1, 0))
    _test_invariant(locus.to_position, 31, locus.to_coordinate, (2, 0))
    _test_invariant(locus.to_position, 33, locus.to_coordinate, (4, 0))
    _test_invariant(locus.to_position, 34, locus.to_coordinate, (5, 0))
    _test_invariant(locus.to_position, 35, locus.to_coordinate, (5, 1))


def test_Locus_inverted():
    locus = Locus(_exons[2], True)

    _test_invariant(locus.to_position, 35, locus.to_coordinate, (1, -1))
    _test_invariant(locus.to_position, 34, locus.to_coordinate, (1, 0))
    _test_invariant(locus.to_position, 33, locus.to_coordinate, (2, 0))
    _test_invariant(locus.to_position, 31, locus.to_coordinate, (4, 0))
    _test_invariant(locus.to_position, 30, locus.to_coordinate, (5, 0))
    _test_invariant(locus.to_position, 29, locus.to_coordinate, (5, 1))


def test_MultiLocus():
    multi_locus = MultiLocus(_exons)

    _test_invariant(
        multi_locus.to_position, 6, multi_locus.to_coordinate, (2, 0))

    _test_invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (10, -1))
    _test_invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (10, 0))
    _test_invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (11, 0))
    _test_invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (13, 0))
    _test_invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (14, 0))
    _test_invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (14, 1))


def test_MultiLocus_adjacent_exons():
    multi_locus = MultiLocus(_adjacent_exons)

    _test_invariant(
        multi_locus.to_position, 2, multi_locus.to_coordinate, (2, 0))
    _test_invariant(
        multi_locus.to_position, 3, multi_locus.to_coordinate, (3, 0))

def test_MultiLocus_offsets_odd():
    multi_locus = MultiLocus([(1, 3), (6, 8)])

    _test_invariant(
        multi_locus.to_position, 3, multi_locus.to_coordinate, (2, 1))
    _test_invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (2, 2))
    _test_invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (3, -1))


def test_MultiLocus_offsets_odd_inverted():
    multi_locus = MultiLocus([(1, 3), (6, 8)], True)

    _test_invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (2, 1))
    _test_invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (2, 2))
    _test_invariant(
        multi_locus.to_position, 3, multi_locus.to_coordinate, (3, -1))


def test_MultiLocus_offsets_even():
    multi_locus = MultiLocus([(1, 3), (7, 9)])

    _test_invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (2, 2))
    _test_invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (3, -2))


def test_MultiLocus_offsets_even_inverted():
    multi_locus = MultiLocus([(1, 3), (7, 9)], True)

    _test_invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (2, 2))
    _test_invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (3, -2))


def test_MultiLocus_inverted():
    multi_locus = MultiLocus(_exons, True)

    _test_invariant(
        multi_locus.to_position, 70, multi_locus.to_coordinate, (2, 0))

    _test_invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (9, -1))
    _test_invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (9, 0))
    _test_invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (10, 0))
    _test_invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (12, 0))
    _test_invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (13, 0))
    _test_invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (13, 1))


def test_MultiLocus_negated():
    multi_locus = MultiLocus(_exons, False, True)

    _test_invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (-10, 1))
    _test_invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (-10, 0))
    _test_invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (-11, 0))
    _test_invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (-13, 0))
    _test_invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (-14, 0))
    _test_invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (-14, -1))


def test_MultiLocus_inverted_negated():
    multi_locus = MultiLocus(_exons, True, True)

    _test_invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (-9, 1))
    _test_invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (-9, 0))
    _test_invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (-10, 0))
    _test_invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (-12, 0))
    _test_invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (-13, 0))
    _test_invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (-13, -1))


def test_Crossmap_genomic():
    crossmap = Crossmap()

    _test_invariant(
        crossmap.coordinate_to_genomic, 0, crossmap.genomic_to_coordinate, 1)
    _test_invariant(
        crossmap.coordinate_to_genomic, 98, crossmap.genomic_to_coordinate, 99)


def test_Crossmap_coding():
    crossmap = Crossmap(_exons, _cds)

    _test_invariant(
        crossmap.coordinate_to_coding, 31,
        crossmap.coding_to_coordinate, (-1, 0, 0))
    _test_invariant(
        crossmap.coordinate_to_coding, 32,
        crossmap.coding_to_coordinate, (1, 0, 1))
    _test_invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (1, 0, 2))


def test_Crossmap_coding_no_utr5():
    crossmap = Crossmap([(10, 20)], (10, 15))

    _test_invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (1, -1, 0))
    _test_invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (1, 0, 1))
    _test_invariant(
        crossmap.coordinate_to_coding, 11,
        crossmap.coding_to_coordinate, (2, 0, 1))

def test_Crossmap_coding_small_utr5():
    crossmap = Crossmap([(10, 20)], (11, 15))

    _test_invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (-1, -1, 0))
    _test_invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (-1, 0, 0))
    _test_invariant(
        crossmap.coordinate_to_coding, 11,
        crossmap.coding_to_coordinate, (1, 0, 1))
    _test_invariant(
        crossmap.coordinate_to_coding, 12,
        crossmap.coding_to_coordinate, (2, 0, 1))

def test_Crossmap_coding_no_utr3():
    crossmap = Crossmap([(10, 20)], (15, 20))

    _test_invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (5, 0, 1))
    _test_invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (5, 1, 2))
    _test_invariant(
        crossmap.coordinate_to_coding, 21,
        crossmap.coding_to_coordinate, (5, 2, 2))

def test_Crossmap_coding_small_utr3():
    crossmap = Crossmap([(10, 20)], (15, 19))

    _test_invariant(
        crossmap.coordinate_to_coding, 18,
        crossmap.coding_to_coordinate, (4, 0, 1))
    _test_invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (1, 0, 2))
    _test_invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (1, 1, 2))
    _test_invariant(
        crossmap.coordinate_to_coding, 21,
        crossmap.coding_to_coordinate, (1, 2, 2))

def test_Crossmap_coding_inverted():
    crossmap = Crossmap(_exons, _cds, True)

    _test_invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (-1, 0, 0))
    _test_invariant(
        crossmap.coordinate_to_coding, 42,
        crossmap.coding_to_coordinate, (1, 0, 1))
    _test_invariant(
        crossmap.coordinate_to_coding, 31,
        crossmap.coding_to_coordinate, (1, 0, 2))
