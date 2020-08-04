from mutalyzer_crossmapper import MultiLocus
from mutalyzer_crossmapper.multi_locus import _offsets

from helper import degenerate_equal, invariant


_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_adjacent_exons = [(1, 3), (3, 5)]


def test_offsets():
    assert _offsets(_exons) == [0, 3, 9, 14, 18, 20]
    assert _offsets(_exons, True) == [0, 2, 4, 8, 13, 19]
    assert _offsets(_adjacent_exons) == [0, 2]
    assert _offsets(_adjacent_exons, True) == [0, 2]


def test_MultiLocus():
    multi_locus = MultiLocus(_exons)

    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (1, -1, -1))
    invariant(
        multi_locus.to_position, 6, multi_locus.to_coordinate, (2, 0, 0))
    invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (10, -1, 0))
    invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (10, 0, 0))
    invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (11, 0, 0))
    invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (13, 0, 0))
    invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (14, 0, 0))
    invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (14, 1, 0))
    invariant(
        multi_locus.to_position, 72, multi_locus.to_coordinate, (22, 1, 1))


def test_MultiLocus_inverted():
    multi_locus = MultiLocus(_exons, True)

    invariant(
        multi_locus.to_position, 70, multi_locus.to_coordinate, (2, 0, 0))

    invariant(
        multi_locus.to_position, 72, multi_locus.to_coordinate, (1, -1, -1))
    invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (9, -1, 0))
    invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (9, 0, 0))
    invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (10, 0, 0))
    invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (12, 0, 0))
    invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (13, 0, 0))
    invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (13, 1, 0))
    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (22, 1, 1))


def test_MultiLocus_negated():
    multi_locus = MultiLocus(_exons, False, True)

    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (-1, 1, 1))
    invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (-10, 1, 0))
    invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (-10, 0, 0))
    invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (-11, 0, 0))
    invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (-13, 0, 0))
    invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (-14, 0, 0))
    invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (-14, -1, 0))
    invariant(
        multi_locus.to_position, 72, multi_locus.to_coordinate, (-22, -1, -1))


def test_MultiLocus_inverted_negated():
    multi_locus = MultiLocus(_exons, True, True)

    invariant(
        multi_locus.to_position, 72, multi_locus.to_coordinate, (-1, 1, 1))
    invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (-9, 1, 0))
    invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (-9, 0, 0))
    invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (-10, 0, 0))
    invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (-12, 0, 0))
    invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (-13, 0, 0))
    invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (-13, -1, 0))
    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (-22, -1, -1))


def test_MultiLocus_adjacent_exons():
    multi_locus = MultiLocus(_adjacent_exons)

    invariant(
        multi_locus.to_position, 0, multi_locus.to_coordinate, (1, -1, -1))
    invariant(
        multi_locus.to_position, 2, multi_locus.to_coordinate, (2, 0, 0))
    invariant(
        multi_locus.to_position, 3, multi_locus.to_coordinate, (3, 0, 0))
    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (4, 1, 1))


def test_MultiLocus_offsets_odd():
    multi_locus = MultiLocus([(1, 3), (6, 8)])

    invariant(
        multi_locus.to_position, 0, multi_locus.to_coordinate, (1, -1, -1))
    invariant(
        multi_locus.to_position, 3, multi_locus.to_coordinate, (2, 1, 0))
    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (2, 2, 0))
    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (3, -1, 0))
    invariant(
        multi_locus.to_position, 8, multi_locus.to_coordinate, (4, 1, 1))


def test_MultiLocus_offsets_odd_inverted():
    multi_locus = MultiLocus([(1, 3), (6, 8)], True)

    invariant(
        multi_locus.to_position, 8, multi_locus.to_coordinate, (1, -1, -1))
    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (2, 1, 0))
    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (2, 2, 0))
    invariant(
        multi_locus.to_position, 3, multi_locus.to_coordinate, (3, -1, 0))
    invariant(
        multi_locus.to_position, 0, multi_locus.to_coordinate, (4, 1, 1))


def test_MultiLocus_offsets_even():
    multi_locus = MultiLocus([(1, 3), (7, 9)])

    invariant(
        multi_locus.to_position, 0, multi_locus.to_coordinate, (1, -1, -1))
    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (2, 2, 0))
    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (3, -2, 0))
    invariant(
        multi_locus.to_position, 9, multi_locus.to_coordinate, (4, 1, 1))


def test_MultiLocus_offsets_even_inverted():
    multi_locus = MultiLocus([(1, 3), (7, 9)], True)

    invariant(
        multi_locus.to_position, 0, multi_locus.to_coordinate, (4, 1, 1))
    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (2, 2, 0))
    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (3, -2, 0))
    invariant(
        multi_locus.to_position, 9, multi_locus.to_coordinate, (1, -1, -1))


def test_MultiLocus_degenerate():
    multi_locus = MultiLocus(_exons)

    degenerate_equal(
        multi_locus.to_coordinate, 4, [(1, -1, -1), (-1, 0, -1)])
    degenerate_equal(
        multi_locus.to_coordinate, 72, [(22, 1, 1), (23, 0, 1)])


def test_MultiLocus_degenerate_return():
    multi_locus = MultiLocus(_exons)

    assert multi_locus.to_position(4, True) == (-1, 0, -1)
    assert multi_locus.to_position(72, True) == (23, 0, 1)


def test_MultiLocus_degenerate_inverted():
    multi_locus = MultiLocus(_exons, True)

    degenerate_equal(
        multi_locus.to_coordinate, 72, [(1, -1, -1), (-1, 0, -1)])
    degenerate_equal(
        multi_locus.to_coordinate, 4, [(22, 1, 1), (23, 0, 1)])


def test_MultiLocus_degenerate_inverted_return():
    multi_locus = MultiLocus(_exons, True)

    assert multi_locus.to_position(72, True) == (-1, 0, -1)
    assert multi_locus.to_position(4, True) == (23, 0, 1)


def test_MultiLocus_degenerate_inverted_negated():
    multi_locus = MultiLocus(_exons, True, True)

    degenerate_equal(
        multi_locus.to_coordinate, 72, [(-1, 1, 1), (1, 0, 1)])
    degenerate_equal(
        multi_locus.to_coordinate, 4, [(-22, -1, -1), (-23, 0, -1)])


def test_MultiLocus_degenerate_inverted_negated_return():
    multi_locus = MultiLocus(_exons, True, True)

    assert multi_locus.to_position(72, True) == (1, 0, 1)
    assert multi_locus.to_position(4, True) == (-23, 0, -1)
