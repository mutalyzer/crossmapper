from mutalyzer_crossmapper import MultiLocus
from mutalyzer_crossmapper.multi_locus import _offsets

from helper import degenerate_equal, invariant

_locations = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]


def test_offsets():
    """Cummulative location lengths."""
    assert _offsets(_locations, 1) == [0, 3, 9, 14, 18, 20]


def test_offsets_inverted():
    """Cummulative location lengths for inverted list of locations."""
    assert _offsets(_locations, -1) == [0, 2, 4, 8, 13, 19]


def test_offsets_adjacent():
    """Cummulative location lengths for adjacent locations."""
    assert _offsets([(1, 3), (3, 5)], 1) == [0, 2]


def test_offsets_adjacent_inverted():
    """Cummulative location lengths for inverted list of adjacent locations."""
    assert _offsets([(1, 3), (3, 5)], -1) == [0, 2]


def test_MultiLocus():
    """Forward oriented MultiLocus."""
    multi_locus = MultiLocus(_locations)

    # Boundary between upstream and the first locus.
    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (0, -1, -1))
    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (0, 0, 0))

    # Internal locus.
    invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (9, -1, 0))
    invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (9, 0, 0))
    invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (10, 0, 0))
    invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (12, 0, 0))
    invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (13, 0, 0))
    invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (13, 1, 0))

    # Boundary between the last locus and downstream.
    invariant(
        multi_locus.to_position, 71, multi_locus.to_coordinate, (21, 0, 0))
    invariant(
        multi_locus.to_position, 72, multi_locus.to_coordinate, (21, 1, 1))


def test_MultiLocus_inverted():
    """Reverse oriented MultiLocus."""
    multi_locus = MultiLocus(_locations, True)

    # Boundary between upstream and the first locus.
    invariant(
        multi_locus.to_position, 72, multi_locus.to_coordinate, (0, -1, -1))
    invariant(
        multi_locus.to_position, 71, multi_locus.to_coordinate, (0, 0, 0))

    # Internal locus.
    invariant(
        multi_locus.to_position, 35, multi_locus.to_coordinate, (8, -1, 0))
    invariant(
        multi_locus.to_position, 34, multi_locus.to_coordinate, (8, 0, 0))
    invariant(
        multi_locus.to_position, 33, multi_locus.to_coordinate, (9, 0, 0))
    invariant(
        multi_locus.to_position, 31, multi_locus.to_coordinate, (11, 0, 0))
    invariant(
        multi_locus.to_position, 30, multi_locus.to_coordinate, (12, 0, 0))
    invariant(
        multi_locus.to_position, 29, multi_locus.to_coordinate, (12, 1, 0))

    # Boundary between the last locus and downstream.
    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (21, 0, 0))
    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (21, 1, 1))


def test_MultiLocus_adjacent_loci():
    """Positions are continuous when loci are adjacent."""
    multi_locus = MultiLocus([(1, 3), (3, 5)])

    invariant(
        multi_locus.to_position, 2, multi_locus.to_coordinate, (1, 0, 0))
    invariant(
        multi_locus.to_position, 3, multi_locus.to_coordinate, (2, 0, 0))


def test_MultiLocus_adjacent_loci_inverted():
    """Positions are continuous when loci are adjacent."""
    multi_locus = MultiLocus([(1, 3), (3, 5)], True)

    invariant(
        multi_locus.to_position, 3, multi_locus.to_coordinate, (1, 0, 0))
    invariant(
        multi_locus.to_position, 2, multi_locus.to_coordinate, (2, 0, 0))


def test_MultiLocus_offsets_odd():
    """Offets exacly between two loci are assigned to the upstream locus."""
    multi_locus = MultiLocus([(1, 3), (6, 8)])

    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (1, 2, 0))
    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (2, -1, 0))


def test_MultiLocus_offsets_odd_inverted():
    """Offets exacly between two loci are assigned to the upstream locus."""
    multi_locus = MultiLocus([(1, 3), (6, 8)], True)

    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (1, 2, 0))
    invariant(
        multi_locus.to_position, 3, multi_locus.to_coordinate, (2, -1, 0))


def test_MultiLocus_offsets_even():
    """Offsets are assigned to the nearest locus."""
    multi_locus = MultiLocus([(1, 3), (7, 9)])

    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (1, 2, 0))
    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (2, -2, 0))


def test_MultiLocus_offsets_even_inverted():
    """Offsets are assigned to the nearest locus."""
    multi_locus = MultiLocus([(1, 3), (7, 9)], True)

    invariant(
        multi_locus.to_position, 5, multi_locus.to_coordinate, (1, 2, 0))
    invariant(
        multi_locus.to_position, 4, multi_locus.to_coordinate, (2, -2, 0))


def test_MultiLocus_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    multi_locus = MultiLocus(_locations)

    degenerate_equal(
        multi_locus.to_coordinate, 4, [(0, -1, -1), (-1, 0, -1)])
    degenerate_equal(
        multi_locus.to_coordinate, 72, [(21, 1, 1), (22, 0, 1)])


def test_MultiLocus_inverted_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    multi_locus = MultiLocus(_locations, True)

    degenerate_equal(
        multi_locus.to_coordinate, 72, [(0, -1, -1), (-1, 0, -1)])
    degenerate_equal(
        multi_locus.to_coordinate, 4, [(21, 1, 1), (22, 0, 1)])
