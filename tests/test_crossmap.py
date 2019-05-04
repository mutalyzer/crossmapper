from crossmapper.crossmapper import (
    Crossmap, Locus, MultiLocus, _offsets)


_exons = [(5, 8), (14, 20), (30, 35), (40, 46), (70, 72)]
_cds = (32, 43)

_offsets_exon = [0, 3, 9, 14, 20]
_offsets_exon_inverted = [0, 2, 8, 13, 19]


def _test_invariant(f, x, f_i, y):
    assert f(x) == y
    assert f_i(y) == x


def test_offsets():
    assert _offsets(_exons) == _offsets_exon
    assert _offsets(_exons, True) == _offsets_exon_inverted


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

    #raise ValueError(crossmap._coding[0]._offsets)
    assert crossmap.coordinate_to_coding(31) == (-1, 0, 0)
    assert crossmap.coordinate_to_coding(32) == (1, 0, 1)
    assert crossmap.coordinate_to_coding(43) == (1, 0, 2)

    _test_invariant(
        crossmap.coordinate_to_coding, 31,
        crossmap.coding_to_coordinate, (-1, 0, 0))
    _test_invariant(
        crossmap.coordinate_to_coding, 32,
        crossmap.coding_to_coordinate, (1, 0, 1))
    _test_invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (1, 0, 2))
