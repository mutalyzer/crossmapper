from crossmapper.crossmapper import (
    _coordinate_to_exon,
    _coordinate_to_locus,
    _exon_to_coordinate,
    _locus_to_coordinate,
    _offsets,
    )


_exons = [(5, 8), (14, 20), (30, 35), (40, 46), (70, 72)]
_cds = (32, 43)

_offsets_exon = [0, 3, 9, 14, 20]
_offsets_exon_inverted = [19, 13, 8, 2, 0]


def _test_invariant(f, x, f_i, y, args):
    assert f(x, *args) == y
    assert f_i(y, *args) == x


def test_offsets():
    assert _offsets(_exons) == _offsets_exon
    assert _offsets(_exons, True) == _offsets_exon_inverted


def test_locus():
    locus = _exons[2]

    _test_invariant(
        _coordinate_to_locus, 29, _locus_to_coordinate, (1, -1), [locus])
    _test_invariant(
        _coordinate_to_locus, 30, _locus_to_coordinate, (1, 0), [locus])
    _test_invariant(
        _coordinate_to_locus, 31, _locus_to_coordinate, (2, 0), [locus])
    _test_invariant(
        _coordinate_to_locus, 33, _locus_to_coordinate, (4, 0), [locus])
    _test_invariant(
        _coordinate_to_locus, 34, _locus_to_coordinate, (5, 0), [locus])
    _test_invariant(
        _coordinate_to_locus, 35, _locus_to_coordinate, (5, 1), [locus])


def test_locus_inverted():
    locus = _exons[2]

    _test_invariant(
        _coordinate_to_locus, 35, _locus_to_coordinate, (1, -1), [locus, True])
    _test_invariant(
        _coordinate_to_locus, 34, _locus_to_coordinate, (1, 0), [locus, True])
    _test_invariant(
        _coordinate_to_locus, 33, _locus_to_coordinate, (2, 0), [locus, True])
    _test_invariant(
        _coordinate_to_locus, 31, _locus_to_coordinate, (4, 0), [locus, True])
    _test_invariant(
        _coordinate_to_locus, 30, _locus_to_coordinate, (5, 0), [locus, True])
    _test_invariant(
        _coordinate_to_locus, 29, _locus_to_coordinate, (5, 1), [locus, True])


def test_exon():
    exon = 2

    _test_invariant(
        _coordinate_to_exon, 29, _exon_to_coordinate, (10, -1),
        [_exons[exon], _offsets_exon[exon]])
    _test_invariant(
        _coordinate_to_exon, 30, _exon_to_coordinate, (10, 0),
        [_exons[exon], _offsets_exon[exon]])
    _test_invariant(
        _coordinate_to_exon, 31, _exon_to_coordinate, (11, 0),
        [_exons[exon], _offsets_exon[exon]])
    _test_invariant(
        _coordinate_to_exon, 33, _exon_to_coordinate, (13, 0),
        [_exons[exon], _offsets_exon[exon]])
    _test_invariant(
        _coordinate_to_exon, 34, _exon_to_coordinate, (14, 0),
        [_exons[exon], _offsets_exon[exon]])
    _test_invariant(
        _coordinate_to_exon, 35, _exon_to_coordinate, (14, 1),
        [_exons[exon], _offsets_exon[exon]])


def test_exon_inverted():
    exon = 2

    _test_invariant(
        _coordinate_to_exon, 35, _exon_to_coordinate, (9, -1),
        [_exons[exon], _offsets_exon_inverted[exon], True])
    _test_invariant(
        _coordinate_to_exon, 34, _exon_to_coordinate, (9, 0),
        [_exons[exon], _offsets_exon_inverted[exon], True])
    _test_invariant(
        _coordinate_to_exon, 33, _exon_to_coordinate, (10, 0),
        [_exons[exon], _offsets_exon_inverted[exon], True])
    _test_invariant(
        _coordinate_to_exon, 31, _exon_to_coordinate, (12, 0),
        [_exons[exon], _offsets_exon_inverted[exon], True])
    _test_invariant(
        _coordinate_to_exon, 30, _exon_to_coordinate, (13, 0),
        [_exons[exon], _offsets_exon_inverted[exon], True])
    _test_invariant(
        _coordinate_to_exon, 29, _exon_to_coordinate, (13, 1),
        [_exons[exon], _offsets_exon_inverted[exon], True])
