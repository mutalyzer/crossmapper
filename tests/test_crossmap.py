from crossmapper.crossmapper import (
    #_coding_offsets,
    _coordinate_to_locus,
    _coordinate_to_noncoding,
    _locus_to_coordinate,
    _noncoding_offsets,
    _noncoding_to_coordinate,
    )


_exons = [(5, 8), (14, 20), (30, 35), (40, 46), (70, 72)]
_cds = (32, 43)

_offsets_noncoding = [0, 3, 9, 14, 20]
_offsets_noncoding_inverted = [19, 13, 8, 2, 0]
#_offsets_coding = [(0, 11), (0, 8), (0, 2), (3, 0), (5, 0)]


def _test_invariant(f, x, f_i, y, args):
    assert f(x, *args) == y
    assert f_i(y, *args) == x


def test_noncoding_offsets():
    assert _noncoding_offsets(_exons) == _offsets_noncoding
    assert _noncoding_offsets(_exons, True) == _offsets_noncoding_inverted


#def test_coding_offsets():
#    assert _coding_offsets(_exons, _cds) == _offsets_coding


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


def test_noncoding():
    exon = 2

    _test_invariant(
        _coordinate_to_noncoding, 29, _noncoding_to_coordinate, (10, -1),
        [_exons[exon], _offsets_noncoding[exon]])
    _test_invariant(
        _coordinate_to_noncoding, 30, _noncoding_to_coordinate, (10, 0),
        [_exons[exon], _offsets_noncoding[exon]])
    _test_invariant(
        _coordinate_to_noncoding, 31, _noncoding_to_coordinate, (11, 0),
        [_exons[exon], _offsets_noncoding[exon]])
    _test_invariant(
        _coordinate_to_noncoding, 33, _noncoding_to_coordinate, (13, 0),
        [_exons[exon], _offsets_noncoding[exon]])
    _test_invariant(
        _coordinate_to_noncoding, 34, _noncoding_to_coordinate, (14, 0),
        [_exons[exon], _offsets_noncoding[exon]])
    _test_invariant(
        _coordinate_to_noncoding, 35, _noncoding_to_coordinate, (14, 1),
        [_exons[exon], _offsets_noncoding[exon]])


def test_noncoding_inverted():
    exon = 2

    _test_invariant(
        _coordinate_to_noncoding, 35, _noncoding_to_coordinate, (9, -1),
        [_exons[exon], _offsets_noncoding_inverted[exon], True])
    _test_invariant(
        _coordinate_to_noncoding, 34, _noncoding_to_coordinate, (9, 0),
        [_exons[exon], _offsets_noncoding_inverted[exon], True])
    _test_invariant(
        _coordinate_to_noncoding, 33, _noncoding_to_coordinate, (10, 0),
        [_exons[exon], _offsets_noncoding_inverted[exon], True])
    _test_invariant(
        _coordinate_to_noncoding, 31, _noncoding_to_coordinate, (12, 0),
        [_exons[exon], _offsets_noncoding_inverted[exon], True])
    _test_invariant(
        _coordinate_to_noncoding, 30, _noncoding_to_coordinate, (13, 0),
        [_exons[exon], _offsets_noncoding_inverted[exon], True])
    _test_invariant(
        _coordinate_to_noncoding, 29, _noncoding_to_coordinate, (13, 1),
        [_exons[exon], _offsets_noncoding_inverted[exon], True])
