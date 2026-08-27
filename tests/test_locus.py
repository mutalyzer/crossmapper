from mutalyzer_crossmapper.locus import Locus, Point, Coord

from helper import invariant
import pytest


def test_invalid_Locus_initialization():
    """Test Locus initialization."""
    with pytest.raises(ValueError) as e:
        Locus((10, 5))
    assert str(e.value) == f"Locus start 10 must be smaller than locus end 5."
    with pytest.raises(ValueError) as e:
        Locus((10, 20, 30))
    assert str(e.value) == f"Locus must be a tuple of two values."
    with pytest.raises(ValueError) as e:
        Locus((10, -5))
    assert str(e.value) == f"Value must be non-negative."
    with pytest.raises(ValueError) as e:
        Locus((10, 20.5))
    assert str(e.value) == f"Value must be an integer."
    with pytest.raises(ValueError) as e:
        Locus((10, None))
    assert str(e.value) == f"Value must be an integer."
    with pytest.raises(ValueError) as e:
        Locus(("10", "20"))
    with pytest.raises(ValueError) as e:
        Locus((10, 10))
    assert str(e.value) == f"Locus start 10 must be smaller than locus end 10."

    #Inverted Locus initialization
    with pytest.raises(ValueError) as e:
        Locus((10, 5), inverted=True)
    assert str(e.value) == f"Locus start 10 must be smaller than locus end 5."
    with pytest.raises(ValueError) as e:
        Locus((10, 20, 30), inverted=True)
    assert str(e.value) == f"Locus must be a tuple of two values."
    with pytest.raises(ValueError) as e:
        Locus((10, -5), inverted=True)
    assert str(e.value) == f"Value must be non-negative."
    with pytest.raises(ValueError) as e:
        Locus((10, 20.5), inverted=True)
    assert str(e.value) == f"Value must be an integer."
    with pytest.raises(ValueError) as e:
        Locus((10, None), inverted=True)
    assert str(e.value) == f"Value must be an integer."
    with pytest.raises(ValueError) as e:
        Locus(("10", "20"), inverted=True)
    assert str(e.value) == f"Value must be an integer."
    with pytest.raises(ValueError) as e:
        Locus((10, 10), inverted=True)
    assert str(e.value) == f"Locus start 10 must be smaller than locus end 10."


def test_invalid_Coord_initialization():
    """Test Coord initialization."""
    with pytest.raises(ValueError) as e:
        Coord(-1)
    assert str(e.value) == f"Value must be non-negative."
    with pytest.raises(ValueError) as e:
        Coord(3.5)
    assert str(e.value) == f"Value must be an integer."
    with pytest.raises(ValueError) as e:
        Coord("10")
    assert str(e.value) == f"Value must be an integer."
    with pytest.raises(ValueError) as e:
        Coord(None)
    assert str(e.value) == f"Value must be an integer."
    with pytest.raises(ValueError) as e:
        Coord([10])
    assert str(e.value) == f"Value must be an integer."


def test_invalid_Locus_point():
    """Forward orientent Locus with invalid point."""
    locus = Locus((30, 35))
    with pytest.raises(ValueError) as e:
        locus.to_coordinate(Point(position=-5, offset=0))
    assert str(e.value) == f"Value must be non-negative."

    with pytest.raises(IndexError) as e:
        locus.to_coordinate(Point(position=5, offset=0))
    assert str(e.value) == f"Position 5 exceeds locus length."
    with pytest.raises(IndexError) as e:
        locus.to_coordinate(Point(position=0, offset=2))
    assert str(e.value) == f"Offset 2 should be at a locus end."
    with pytest.raises(IndexError) as e:
        locus.to_coordinate(Point(position=4, offset=-2))
    assert str(e.value) == f"Offset -2 should be at a locus start."
    with pytest.raises(ValueError) as e:
        locus.to_coordinate(Point(position=2, offset=1))
    assert str(e.value) == f"Position 2 is not at a locus boundary."


def test_invalid_Locus_inverted_point():
    """Reverse orientent Locus with invalid point."""
    locus = Locus((30, 35), True)
    with pytest.raises(ValueError) as e:
        locus.to_coordinate(Point(position=-5, offset=0))
    assert str(e.value) == f"Value must be non-negative."
    with pytest.raises(IndexError) as e:
        locus.to_coordinate(Point(position=5, offset=0))
    assert str(e.value) == f"Position 5 exceeds locus length."
    with pytest.raises(IndexError) as e:
        locus.to_coordinate(Point(position=0, offset=2))
    assert str(e.value) == f"Offset 2 should be at a locus end."
    with pytest.raises(IndexError) as e:
        locus.to_coordinate(Point(position=4, offset=-2))
    assert str(e.value) == f"Offset -2 should be at a locus start."
    with pytest.raises(ValueError) as e:
        locus.to_coordinate(Point(position=2, offset=1))
    assert str(e.value) == f"Position 2 is not at a locus boundary."


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
