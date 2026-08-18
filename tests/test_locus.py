from mutalyzer_crossmapper.locus import Locus, Point, Coord

from helper import invariant
import pytest


def test_invalid_Locus_initialization():
    """Test Locus initialization."""
    with pytest.raises(ValueError):
        Locus((10, 5))
    with pytest.raises(ValueError):
        Locus((10, 20, 30))
    with pytest.raises(ValueError):
        Locus((10, -5))
    with pytest.raises(ValueError):
        Locus((10, 20.5))
    with pytest.raises(ValueError):
        Locus((10.5, None))
    with pytest.raises(ValueError):
        Locus(("10", "20"))
    with pytest.raises(ValueError):
        Locus((10, 10))

    #Inverted Locus initialization
    with pytest.raises(ValueError):
        Locus((10, 5), inverted=True)
    with pytest.raises(ValueError):
        Locus((10, 20, 30), inverted=True)
    with pytest.raises(ValueError):
        Locus((10, -5), inverted=True)
    with pytest.raises(ValueError):
        Locus((10, 20.5), inverted=True)
    with pytest.raises(ValueError):
        Locus((10.5, None), inverted=True)
    with pytest.raises(ValueError):
        Locus(("10", "20"), inverted=True)
    with pytest.raises(ValueError):
        Locus((10, 10), inverted=True)


def test_invalid_Coord_initialization():
    """Test Coord initialization."""
    with pytest.raises(ValueError):
        Coord(-1)
    with pytest.raises(ValueError):
        Coord(3.5)
    with pytest.raises(ValueError):
        Coord("10")
    with pytest.raises(ValueError):
        Coord(None)
    with pytest.raises(ValueError):
        Coord([10])


def test_invalid_Locus_point():
    """Forward orientent Locus with invalid point."""
    locus = Locus((30, 35))
    with pytest.raises(ValueError):
        locus.to_coordinate(Point(position=-5, offset=0))
    with pytest.raises(IndexError):
        locus.to_coordinate(Point(position=5, offset=0))
    with pytest.raises(IndexError):
        locus.to_coordinate(Point(position=0, offset=2))
    with pytest.raises(IndexError):
        locus.to_coordinate(Point(position=4, offset=-2))
    with pytest.raises(ValueError):
        locus.to_coordinate(Point(position=2, offset=1))


def test_invalid_Locus_inverted_point():
    """Reverse orientent Locus with invalid point."""
    locus = Locus((30, 35), True)
    with pytest.raises(ValueError):
        locus.to_coordinate(Point(position=-5, offset=0))
    with pytest.raises(IndexError):
        locus.to_coordinate(Point(position=5, offset=0))
    with pytest.raises(IndexError):
        locus.to_coordinate(Point(position=0, offset=2))
    with pytest.raises(IndexError):
        locus.to_coordinate(Point(position=4, offset=-2))
    with pytest.raises(ValueError):
        locus.to_coordinate(Point(position=2, offset=1))


def test_Locus():
    """Forward orientent Locus."""
    locus = Locus((30, 35))

    invariant(locus.to_position, Coord(29), locus.to_coordinate, Point(position=0, offset=-1))
    invariant(locus.to_position, Coord(30), locus.to_coordinate, Point(position=0, offset=0))
    invariant(locus.to_position, Coord(31), locus.to_coordinate, Point(position=1, offset=0))
    invariant(locus.to_position, Coord(33), locus.to_coordinate, Point(position=3, offset=0))
    invariant(locus.to_position, Coord(34), locus.to_coordinate, Point(position=4, offset=0))
    invariant(locus.to_position, Coord(35), locus.to_coordinate, Point(position=4, offset=1))


def test_Locus_inverted():
    """Reverse orientent Locus."""
    locus = Locus((30, 35), True)

    invariant(locus.to_position, Coord(35), locus.to_coordinate, Point(position=0, offset=-1))
    invariant(locus.to_position, Coord(34), locus.to_coordinate, Point(position=0, offset=0))
    invariant(locus.to_position, Coord(33), locus.to_coordinate, Point(position=1, offset=0))
    invariant(locus.to_position, Coord(31), locus.to_coordinate, Point(position=3, offset=0))
    invariant(locus.to_position, Coord(30), locus.to_coordinate, Point(position=4, offset=0))
    invariant(locus.to_position, Coord(29), locus.to_coordinate, Point(position=4, offset=1))
