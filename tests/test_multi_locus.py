from mutalyzer_crossmapper import multi_locus
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
    with pytest.raises(ValueError) as e:
        MultiLocus(([(10, 5), (20, 25)]))
    assert str(e.value) == 'Locus start 10 must be smaller than locus end 5.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, 20, 30), (40, 50)])
    assert str(e.value) == 'Locus must be a tuple of two values.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, -5), (20, 25)])
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, 20.5), (30, 40)])
    assert str(e.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10.5, None), (20, 30)])
    assert str(e.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as e:
        MultiLocus([('10', '20'), (30, 40)])
    assert str(e.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, 20), (15, 25)])
    assert str(e.value) == 'Locus (15, 25) and locus (10, 20) are overlapping.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, 12), (15, 25)], 25)
    assert str(e.value) == 'Value 25 must be within the bounds of the reference length 25.'

    # Inverted MultiLocus initialization
    with pytest.raises(ValueError) as e:
        MultiLocus(([(10, 5), (20, 25)]), inverted=True)
    assert str(e.value) == 'Locus start 10 must be smaller than locus end 5.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, 20, 30), (40, 50)], inverted=True)
    assert str(e.value) == 'Locus must be a tuple of two values.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, -5), (20, 25)], inverted=True)
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, 20.5), (30, 40)], inverted=True)
    assert str(e.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10.5, None), (20, 30)], inverted=True)
    assert str(e.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as e:
        MultiLocus([('10', '20'), (30, 40)], inverted=True)
    assert str(e.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, 20), (15, 25)], 25, inverted=True)
    assert str(e.value) == 'Locus (15, 25) and locus (10, 20) are overlapping.'
    with pytest.raises(ValueError) as e:
        MultiLocus([(10, 12), (15, 25)], 25, inverted=True)
    assert str(e.value) == 'Value 25 must be within the bounds of the reference length 25.'


def test_MultiLocus_invalid_coordinate():
    """Forward orientent MultiLocus with invalid coordinate."""
    multi_locus = MultiLocus([(30, 35), (40, 45)])
    with pytest.raises(ValueError) as e:
        multi_locus.to_position(Coord(-1))
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_position(Coord(46.7))
    assert str(e.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_position(Coord('31'))
    assert str(e.value) == 'Value must be an integer.'


def test_invalid_Point_initialization():
    """Test Point initialization."""
    with pytest.raises(ValueError) as e:
        Point(position=-1, offset=0, region='')
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        Point(position=0, offset=0, region='*')
    assert str(e.value) == 'Region * is not valid. Must be "", "u", or "d".'
    with pytest.raises(ValueError) as e:
        Point(position=0, offset=0, region=None)
    assert str(e.value) == 'Region None is not valid. Must be "", "u", or "d".'
    with pytest.raises(ValueError) as e:
        Point(position='11', offset=0, region='u')
    assert str(e.value) == 'Value must be an integer.'


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
    with pytest.raises(ValueError) as e:
        multi_locus.to_position(Coord(74))
    with pytest.raises(ValueError) as e:
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
    with pytest.raises(ValueError) as e:
        multi_locus.to_position(Coord(74))
    with pytest.raises(ValueError) as e:
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
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=1, offset=-1, region='u'))
    assert str(e.value) == 'Position 1 is not at upstream boundary.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=-1, offset=-1, region='u'))
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=20, offset=-1, region='u'))
    assert str(e.value) == 'Position 20 is not at upstream boundary.'


def test_upstream_invalid_position_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=1, offset=-1, region='u'))
    assert str(e.value) == 'Position 1 is not at upstream boundary.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=-1, offset=-1, region='u'))
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=20, offset=-1, region='u'))
    assert str(e.value) == 'Position 20 is not at upstream boundary.'


def test_upstream_invalid_offset():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=1, region='u'))
    assert str(e.value) == 'Offset 1 at upstream boundary should be negative.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=-6, region='u'))
    assert str(e.value) == 'Offset -6 exceeds upstream region.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=0, region='u'))
    assert str(e.value) == 'Offset 0 at upstream boundary should be negative.'


def test_upstream_invalid_offset_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=1, region='u'))
    assert str(e.value) == 'Offset 1 at upstream boundary should be negative.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=-6, region='u'))
    assert str(e.value) == 'Offset -6 exceeds upstream region.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=0, region='u'))
    assert str(e.value) == 'Offset 0 at upstream boundary should be negative.'


def test_transcribed_invalid_position():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=-1, offset=0, region=''))
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=7, offset=1, region=''))
    assert str(e.value) == 'Position 7 is not at a locus boundary.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=20, offset=0, region=''))
    assert str(e.value) == 'Position 20 exceeds multi locus length.'


def test_transcribed_invalid_position_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=-1, offset=0, region=''))
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=7, offset=-1, region=''))
    assert str(e.value) == 'Position 7 is not at a locus boundary.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=20, offset=0, region=''))
    assert str(e.value) == 'Position 20 exceeds multi locus length.'


def test_transcribed_invalid_offset():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=-5, region=''))
    assert str(e.value) == 'Offset -5 at the first exon should be in the upstream region.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=1, region=''))
    assert str(e.value) == 'Offset 1 should be at a locus end.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=10, region=''))
    assert str(e.value) == 'Offset 10 exceeds intron length.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=4, offset=6, region=''))
    assert str(e.value) == 'Offset 6 exceeds intron length.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=4, offset=-1, region=''))
    assert str(e.value) == 'Offset -1 should be at a locus start.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=5, offset=2, region=''))
    assert str(e.value) == 'Offset 2 should be at a locus end.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=5, offset=-6, region=''))
    assert str(e.value) == 'Offset -6 exceeds intron length.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=1, region=''))
    assert str(e.value) == 'Offset 1 at the last exon should be in the downstream region.'


def test_transcribed_invalid_offset_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=5, region=''))
    assert str(e.value) == 'Offset 5 at the first exon on the reverse complement should be in the downstream region.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=-1, region=''))
    assert str(e.value) == 'Offset -1 should be at a locus start.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=-10, region=''))
    assert str(e.value) == 'Offset -10 exceeds intron length.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=5, offset=-6, region=''))
    assert str(e.value) == 'Offset -6 exceeds intron length.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=5, offset=1, region=''))
    assert str(e.value) == 'Offset 1 should be at a locus end.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=4, offset=6, region=''))
    assert str(e.value) == 'Offset 6 exceeds intron length.'
    with pytest.raises(IndexError) as e:
        multi_locus.to_coordinate(Point(position=4, offset=-1, region=''))
    assert str(e.value) == 'Offset -1 should be at a locus start.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=-5, region=''))
    assert str(e.value) == 'Offset -5 at the first exon on the reverse complement should be in the upstream region.'


def test_downstream_invalid_position():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=1, region='d'))
    assert str(e.value) == 'Position 0 is not at downstream boundary.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=-1, offset=1, region='d'))
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=11, offset=1, region='d'))
    assert str(e.value) == 'Position 11 is not at downstream boundary.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=100, offset=1, region='d'))
    assert str(e.value) == 'Position 100 is not at downstream boundary.'


def test_downstream_invalid_position_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=0, offset=1, region='d'))
    assert str(e.value) == 'Position 0 is not at downstream boundary.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=-1, offset=1, region='d'))
    assert str(e.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=11, offset=1, region='d'))
    assert str(e.value) == 'Position 11 is not at downstream boundary.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=100, offset=1, region='d'))
    assert str(e.value) == 'Position 100 is not at downstream boundary.'


def test_downstream_invalid_offset():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=-1, region='d'))
    assert str(e.value) == 'Offset -1 at downstream boundary should be positive.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=0, region='d'))
    assert str(e.value) == 'Offset 0 at downstream boundary should be positive.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=-5, region='d'))
    assert str(e.value) == 'Offset -5 at downstream boundary should be positive.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=6, region='d'))
    assert str(e.value) == 'Offset 6 exceeds downstream region.'


def test_downstream_invalid_offset_inverted():
    multi_locus = MultiLocus([(5, 10), (15, 20)], length=25, inverted=True)
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=-1, region='d'))
    assert str(e.value) == 'Offset -1 at downstream boundary should be positive.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=0, region='d'))
    assert str(e.value) == 'Offset 0 at downstream boundary should be positive.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=-5, region='d'))
    assert str(e.value) == 'Offset -5 at downstream boundary should be positive.'
    with pytest.raises(ValueError) as e:
        multi_locus.to_coordinate(Point(position=9, offset=6, region='d'))
    assert str(e.value) == 'Offset 6 exceeds downstream region.'

