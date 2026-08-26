from mutalyzer_crossmapper.multi_locus import _offsets, Coord, MultiLocus, Point
from helper import invariant

import pytest

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


## Test MultiLocus model and its point model
def test_invalid_MultiLocus_initialization():
    """Test MultiLocus initialization."""
    with pytest.raises(ValueError):
        MultiLocus(([(10, 5), (20, 25)]))
    with pytest.raises(ValueError):
        MultiLocus([(10, 20, 30), (40, 50)])
    with pytest.raises(ValueError):
        MultiLocus([(10, -5), (20, 25)])
    with pytest.raises(ValueError):
        MultiLocus([(10, 20.5), (30, 40)])
    with pytest.raises(ValueError):
        MultiLocus([(10.5, None), (20, 30)])
    with pytest.raises(ValueError):
        MultiLocus([("10", "20"), (30, 40)])
    with pytest.raises(ValueError):
        MultiLocus([(10, 20), (15, 25)])
    with pytest.raises(ValueError):
        MultiLocus([(10, 12), (15, 25)], 25)

    # Inverted MultiLocus initialization
    with pytest.raises(ValueError):
        MultiLocus(([(10, 5), (20, 25)]), inverted=True)
    with pytest.raises(ValueError):
        MultiLocus([(10, 20, 30), (40, 50)], inverted=True)
    with pytest.raises(ValueError):
        MultiLocus([(10, -5), (20, 25)], inverted=True)
    with pytest.raises(ValueError):
        MultiLocus([(10, 20.5), (30, 40)], inverted=True)
    with pytest.raises(ValueError):
        MultiLocus([(10.5, None), (20, 30)], inverted=True)
    with pytest.raises(ValueError):
        MultiLocus([("10", "20"), (30, 40)], inverted=True)
    with pytest.raises(ValueError):
        MultiLocus([(10, 20), (15, 25)], 25, inverted=True)
    with pytest.raises(ValueError):
        MultiLocus([(10, 12), (15, 25)], 25, inverted=True)


def test_MultiLocus_invalid_coordinate():
    """Forward orientent MultiLocus with invalid coordinate."""
    multi_locus = MultiLocus([(30, 35), (40, 45)])
    with pytest.raises(ValueError):
        multi_locus.to_position(Coord(-1))
    with pytest.raises(ValueError):
        multi_locus.to_position(Coord(46.7))
    with pytest.raises(ValueError):
        multi_locus.to_position(Coord("31"))


def test_invalid_Point_initialization():
    """Test Point initialization."""
    with pytest.raises(ValueError):
        Point(position=-1, offset=0, region='')
    with pytest.raises(ValueError):
        Point(position=0, offset=0, region='*')
    with pytest.raises(ValueError):
        Point(position=0, offset=0, region=None)
    with pytest.raises(ValueError):
        Point(position="11", offset=0, region='u')


def test_MultiLocus():
    """Forward oriented MultiLocus."""
    multi_locus = MultiLocus(_locations)

    # Boundary between upstream and the first locus.
    invariant(
        multi_locus.to_position,
        Coord(4),
        multi_locus.to_coordinate,
        Point(position=0, offset=-1, region='u'),
    )

    invariant(
        multi_locus.to_position,
        Coord(5),
        multi_locus.to_coordinate,
        Point(position=0, offset=0, region=''),
    )

    # Internal locus.
    invariant(
        multi_locus.to_position,
        Coord(29),
        multi_locus.to_coordinate,
        Point(position=9, offset=-1, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(30),
        multi_locus.to_coordinate,
        Point(position=9, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(31),
        multi_locus.to_coordinate,
        Point(position=10, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(33),
        multi_locus.to_coordinate,
        Point(position=12, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(34),
        multi_locus.to_coordinate,
        Point(position=13, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(35),
        multi_locus.to_coordinate,
        Point(position=13, offset=1, region=''),
    )

    # Boundary between the last locus and downstream.
    invariant(
        multi_locus.to_position,
        Coord(71),
        multi_locus.to_coordinate,
        Point(position=21, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(72),
        multi_locus.to_coordinate,
        Point(position=21, offset=1, region='d'),
    )


def test_MultiLocus_inverted():
    """Reverse oriented MultiLocus."""
    multi_locus = MultiLocus(_locations, None, True)

    # Boundary between upstream and the first locus.
    invariant(
        multi_locus.to_position,
        Coord(72),
        multi_locus.to_coordinate,
        Point(position=0, offset=-1, region='u'),
    )
    invariant(
        multi_locus.to_position,
        Coord(71),
        multi_locus.to_coordinate,
        Point(position=0, offset=0, region=''),
    )

    # Internal locus.
    invariant(
        multi_locus.to_position,
        Coord(35),
        multi_locus.to_coordinate,
        Point(position=8, offset=-1, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(34),
        multi_locus.to_coordinate,
        Point(position=8, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(33),
        multi_locus.to_coordinate,
        Point(position=9, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(31),
        multi_locus.to_coordinate,
        Point(position=11, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(30),
        multi_locus.to_coordinate,
        Point(position=12, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(29),
        multi_locus.to_coordinate,
        Point(position=12, offset=1, region=''),
    )

    # Boundary between the last locus and downstream.
    invariant(
        multi_locus.to_position,
        Coord(5),
        multi_locus.to_coordinate,
        Point(position=21, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(4),
        multi_locus.to_coordinate,
        Point(position=21, offset=1, region='d'),
    )


def test_MultiLocus_with_length():
    """Forward oriented MultiLocus."""
    multi_locus = MultiLocus(_locations, length=74)

    # Boundary between the last locus and downstream.
    invariant(
        multi_locus.to_position,
        Coord(71),
        multi_locus.to_coordinate,
        Point(position=21, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(72),
        multi_locus.to_coordinate,
        Point(position=21, offset=1, region='d'),
    )
    invariant(
        multi_locus.to_position,
        Coord(73),
        multi_locus.to_coordinate,
        Point(position=21, offset=2, region='d'),
    )
    # Boundary between the last base and beyond the last base.
    with pytest.raises(ValueError):
        multi_locus.to_position(Coord(74))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=21, offset=3, region='d'))


def test_MultiLocus_inverted_with_length():
    """Inverted MultiLocus with length."""
    multi_locus = MultiLocus(_locations, length=74, inverted=True)

    # Boundary between the first locus and upstream.
    invariant(
        multi_locus.to_position,
        Coord(71),
        multi_locus.to_coordinate,
        Point(position=0, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(72),
        multi_locus.to_coordinate,
        Point(position=0, offset=-1, region='u'),
    )
    # Boundary between the first base beyond the first base.
    invariant(
        multi_locus.to_position,
        Coord(73),
        multi_locus.to_coordinate,
        Point(position=0, offset=-2, region='u'),
    )
    with pytest.raises(ValueError):
        multi_locus.to_position(Coord(74))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=-3, region='u'))


def test_MultiLocus_adjacent_loci():
    """Positions are continuous when loci are adjacent."""
    multi_locus = MultiLocus([(1, 3), (3, 5)])

    invariant(
        multi_locus.to_position,
        Coord(2),
        multi_locus.to_coordinate,
        Point(position=1, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(3),
        multi_locus.to_coordinate,
        Point(position=2, offset=0, region=''),
    )


def test_MultiLocus_adjacent_loci_inverted():
    """Positions are continuous when loci are adjacent."""
    multi_locus = MultiLocus([(1, 3), (3, 5)], None, True)

    invariant(
        multi_locus.to_position,
        Coord(3),
        multi_locus.to_coordinate,
        Point(position=1, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(2),
        multi_locus.to_coordinate,
        Point(position=2, offset=0, region=''),
    )


def test_MultiLocus_offsets_odd():
    """Offets exacly between two loci are assigned to the upstream locus."""
    multi_locus = MultiLocus([(1, 3), (6, 8)])

    invariant(
        multi_locus.to_position,
        Coord(4),
        multi_locus.to_coordinate,
        Point(position=1, offset=2, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(5),
        multi_locus.to_coordinate,
        Point(position=2, offset=-1, region=''),
    )


def test_MultiLocus_offsets_odd_inverted():
    """Offets exacly between two loci are assigned to the upstream locus."""
    multi_locus = MultiLocus([(1, 3), (6, 8)], None, True)
    invariant(
        multi_locus.to_position,
        Coord(4),
        multi_locus.to_coordinate,
        Point(position=1, offset=2, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(3),
        multi_locus.to_coordinate,
        Point(position=2, offset=-1, region=''),
    )


def test_MultiLocus_offsets_even():
    """Offsets are assigned to the nearest locus."""
    multi_locus = MultiLocus([(1, 3), (7, 9)])

    invariant(
        multi_locus.to_position,
        Coord(4),
        multi_locus.to_coordinate,
        Point(position=1, offset=2, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(5),
        multi_locus.to_coordinate,
        Point(position=2, offset=-2, region=''),
    )


def test_MultiLocus_offsets_even_inverted():
    """Offsets are assigned to the nearest locus."""
    multi_locus = MultiLocus([(1, 3), (7, 9)], None, True)

    invariant(
        multi_locus.to_position,
        Coord(5),
        multi_locus.to_coordinate,
        Point(position=1, offset=2, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(4),
        multi_locus.to_coordinate,
        Point(position=2, offset=-2, region=''),
    )


def test_one_base_exon():
    """One base exons."""
    multi_locus = MultiLocus([(1, 2), (4, 5)])
    invariant(
        multi_locus.to_position,
        Coord(0),
        multi_locus.to_coordinate,
        Point(position=0, offset=-1, region='u'),
    )
    invariant(
        multi_locus.to_position,
        Coord(1),
        multi_locus.to_coordinate,
        Point(position=0, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(2),
        multi_locus.to_coordinate,
        Point(position=0, offset=1, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(3),
        multi_locus.to_coordinate,
        Point(position=1, offset=-1, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(4),
        multi_locus.to_coordinate,
        Point(position=1, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(5),
        multi_locus.to_coordinate,
        Point(position=1, offset=1, region='d'),
    )


def test_one_base_exon_inverted():
    """One base exons."""
    multi_locus = MultiLocus([(1, 2), (4, 5)], None, True)
    invariant(
        multi_locus.to_position,
        Coord(0),
        multi_locus.to_coordinate,
        Point(position=1, offset=1, region='d'),
    )
    invariant(
        multi_locus.to_position,
        Coord(1),
        multi_locus.to_coordinate,
        Point(position=1, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(2),
        multi_locus.to_coordinate,
        Point(position=1, offset=-1, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(3),
        multi_locus.to_coordinate,
        Point(position=0, offset=1, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(4),
        multi_locus.to_coordinate,
        Point(position=0, offset=0, region=''),
    )
    invariant(
        multi_locus.to_position,
        Coord(5),
        multi_locus.to_coordinate,
        Point(position=0, offset=-1, region='u'),
    )


def test_upstream_invalid_position():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=1, offset=-1, region='u'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=-1, offset=-1, region='u'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=20, offset=-1, region='u'))


def test_upstream_invalid_position_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=1, offset=-1, region='u'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=-1, offset=-1, region='u'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=20, offset=-1, region='u'))


def test_upstream_invalid_offset():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=1, region='u'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=-6, region='u'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=0, region='u'))


def test_upstream_invalid_offset_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=1, region='u'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=-6, region='u'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=0, region='u'))


def test_transcribed_invalid_position():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=-1, offset=0, region=''))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=7, offset=1, region=''))
    with pytest.raises(IndexError):
        multi_locus.to_coordinate(Point(position=20, offset=0, region=''))


def test_transcribed_invalid_position_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=-1, offset=0, region=''))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=7, offset=-1, region=''))
    with pytest.raises(IndexError):
        multi_locus.to_coordinate(Point(position=20, offset=0, region=''))


def test_transcribed_invalid_offset():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=-5, region=''))
    with pytest.raises(IndexError):
        multi_locus.to_coordinate(Point(position=0, offset=10, region=''))
    with pytest.raises(IndexError):
        multi_locus.to_coordinate(Point(position=4, offset=6, region=''))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=4, offset=-1, region=''))
    with pytest.raises(IndexError):
        multi_locus.to_coordinate(Point(position=5, offset=-6, region=''))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=5, offset=1, region=''))



def test_transcribed_invalid_offset_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=9, offset=5, region=''))
    with pytest.raises(IndexError):
        multi_locus.to_coordinate(Point(position=9, offset=-10, region=''))
    with pytest.raises(IndexError):
        multi_locus.to_coordinate(Point(position=5, offset=-6, region=''))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=5, offset=1, region=''))
    with pytest.raises(IndexError):
        multi_locus.to_coordinate(Point(position=4, offset=6, region=''))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=4, offset=-1, region=''))


def test_downstream_invalid_position():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=1, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=-1, offset=1, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=11, offset=1, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=100, offset=1, region='d'))


def test_downstream_invalid_position_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=0, offset=1, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=-1, offset=1, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=11, offset=1, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=100, offset=1, region='d'))


def test_downstream_invalid_offset():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=9, offset=-1, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=9, offset=0, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=9, offset=-5, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=9, offset=6, region='d'))


def test_downstream_invalid_offset_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=9, offset=-1, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=9, offset=0, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=9, offset=-5, region='d'))
    with pytest.raises(ValueError):
        multi_locus.to_coordinate(Point(position=9, offset=6, region='d'))

