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
        multi_locus.to_position,
        4,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'u'},
    )

    invariant(
        multi_locus.to_position,
        5,
        multi_locus.to_coordinate,
        {'position': 0, 'offset': 0, 'region': ''},
    )

    # Internal locus.
    invariant(
        multi_locus.to_position,
        29,
        multi_locus.to_coordinate,
        {'position': 9, 'offset': -1, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        30,
        multi_locus.to_coordinate,
        {'position': 9, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        31,
        multi_locus.to_coordinate,
        {'position': 10, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        33,
        multi_locus.to_coordinate,
        {'position': 12, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        34,
        multi_locus.to_coordinate,
        {'position': 13, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        35,
        multi_locus.to_coordinate,
        {'position': 13, 'offset': 1, 'region': ''},
    )

    # Boundary between the last locus and downstream.
    invariant(
        multi_locus.to_position,
        71,
        multi_locus.to_coordinate,
        {'position': 21, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        72,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'd'},
    )


def test_MultiLocus_inverted():
    """Reverse oriented MultiLocus."""
    multi_locus = MultiLocus(_locations, True)

    # Boundary between upstream and the first locus.
    invariant(
        multi_locus.to_position,
        72,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'u'},
    )
    invariant(
        multi_locus.to_position,
        71,
        multi_locus.to_coordinate,
        {'position': 0, 'offset': 0, 'region': ''},
    )

    # Internal locus.
    invariant(
        multi_locus.to_position,
        35,
        multi_locus.to_coordinate,
        {'position': 8, 'offset': -1, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        34,
        multi_locus.to_coordinate,
        {'position': 8, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        33,
        multi_locus.to_coordinate,
        {'position': 9, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        31,
        multi_locus.to_coordinate,
        {'position': 11, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        30,
        multi_locus.to_coordinate,
        {'position': 12, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        29,
        multi_locus.to_coordinate,
        {'position': 12, 'offset': 1, 'region': ''},
    )

    # Boundary between the last locus and downstream.
    invariant(
        multi_locus.to_position,
        5,
        multi_locus.to_coordinate,
        {'position': 21, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        4,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'd'},
    )


def test_MultiLocus_adjacent_loci():
    """Positions are continuous when loci are adjacent."""
    multi_locus = MultiLocus([(1, 3), (3, 5)])

    invariant(
        multi_locus.to_position,
        2,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        3,
        multi_locus.to_coordinate,
        {'position': 2, 'offset': 0, 'region': ''},
    )


def test_MultiLocus_adjacent_loci_inverted():
    """Positions are continuous when loci are adjacent."""
    multi_locus = MultiLocus([(1, 3), (3, 5)], True)

    invariant(
        multi_locus.to_position,
        3,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        2,
        multi_locus.to_coordinate,
        {'position': 2, 'offset': 0, 'region': ''},
    )


def test_MultiLocus_offsets_odd():
    """Offets exacly between two loci are assigned to the upstream locus."""
    multi_locus = MultiLocus([(1, 3), (6, 8)])

    invariant(
        multi_locus.to_position,
        4,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 2, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        5,
        multi_locus.to_coordinate,
        {'position': 2, 'offset': -1, 'region': ''},
    )


def test_MultiLocus_offsets_odd_inverted():
    """Offets exacly between two loci are assigned to the upstream locus."""
    multi_locus = MultiLocus([(1, 3), (6, 8)], True)

    invariant(
        multi_locus.to_position,
        4,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 2, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        3,
        multi_locus.to_coordinate,
        {'position': 2, 'offset': -1, 'region': ''},
    )


def test_MultiLocus_offsets_even():
    """Offsets are assigned to the nearest locus."""
    multi_locus = MultiLocus([(1, 3), (7, 9)])

    invariant(
        multi_locus.to_position,
        4,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 2, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        5,
        multi_locus.to_coordinate,
        {'position': 2, 'offset': -2, 'region': ''},
    )


def test_MultiLocus_offsets_even_inverted():
    """Offsets are assigned to the nearest locus."""
    multi_locus = MultiLocus([(1, 3), (7, 9)], True)

    invariant(
        multi_locus.to_position,
        5,
        multi_locus.to_coordinate,
        {'position': 1, 'offset': 2, 'region': ''},
    )
    invariant(
        multi_locus.to_position,
        4,
        multi_locus.to_coordinate,
        {'position': 2, 'offset': -2, 'region': ''},
    )


def test_MultiLocus_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    multi_locus = MultiLocus(_locations)

    degenerate_equal(
        multi_locus.to_coordinate,
        4,
        [
            {'position': 0, 'offset': -1, 'region': 'u'},
            {'position': 1, 'offset': 0, 'region': 'u'},
            {'position': -1, 'offset': 0, 'region': 'u'},
        ],
    )

    degenerate_equal(
        multi_locus.to_coordinate,
        72,
        [
            {'position': 0, 'offset': 1, 'region': 'd'},
            {'position': 1, 'offset': 0, 'region': 'd'},
        ],
    )


def test_MultiLocus_inverted_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    multi_locus = MultiLocus(_locations, True)

    degenerate_equal(
        multi_locus.to_coordinate,
        72,
        [
            {'position': 0, 'offset': 1, 'region': 'u'},
            {'position': -1, 'offset': 0, 'region': 'u'},
            {'position': 1, 'offset': 0, 'region': 'u'},
        ],
    )

    degenerate_equal(
        multi_locus.to_coordinate,
        4,
        [
            {'position': 0, 'offset': -1, 'region': 'd'},
            {'position': 1, 'offset': 0, 'region': 'd'},
        ],
    )
