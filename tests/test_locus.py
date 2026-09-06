import pytest

from helper import invariant
from mutalyzer_crossmapper.locus import Locus, Point


def test_invalid_locus_initialization():
    """Test Locus initialization."""
    with pytest.raises(ValueError) as error:
        Locus((10, 5))
    assert str(error.value) == 'Locus start 10 must be smaller than locus end 5.'
    with pytest.raises(ValueError) as error:
        Locus((10, 20, 30))
    assert str(error.value) == 'Locus must be a tuple of two values.'
    with pytest.raises(ValueError) as error:
        Locus((10, -5))
    assert str(error.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as error:
        Locus((10, 20.5))
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        Locus((10, None))
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        Locus(('10', '20'))
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        Locus((10, True))
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        Locus((10, 10))
    assert str(error.value) == 'Locus start 10 must be smaller than locus end 10.'

    # Inverted Locus initialization
    with pytest.raises(ValueError) as error:
        Locus((10, 5), inverted=True)
    assert str(error.value) == 'Locus start 10 must be smaller than locus end 5.'
    with pytest.raises(ValueError) as error:
        Locus((10, 20, 30), inverted=True)
    assert str(error.value) == 'Locus must be a tuple of two values.'
    with pytest.raises(ValueError) as error:
        Locus((10, -5), inverted=True)
    assert str(error.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as error:
        Locus((10, 20.5), inverted=True)
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        Locus((10, None), inverted=True)
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        Locus(('10', '20'), inverted=True)
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        Locus((10, 10), inverted=True)
    assert str(error.value) == 'Locus start 10 must be smaller than locus end 10.'


def test_invalid_locus_coordinate():
    """Forward orientent Locus with invalid coordinate."""
    locus = Locus((30, 35))

    with pytest.raises(ValueError) as error:
        locus.to_position(-1)
    assert str(error.value) == 'Value must be non-negative.'
    with pytest.raises(ValueError) as error:
        locus.to_position(3.5)
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        locus.to_position('10')
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        locus.to_position(None)
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        locus.to_position([10])
    assert str(error.value) == 'Value must be an integer.'
    with pytest.raises(ValueError) as error:
        locus.to_position(True)
    assert str(error.value) == 'Value must be an integer.'


def test_invalid_locus_point():
    """Forward orientent Locus with invalid point."""
    locus = Locus((30, 35))

    with pytest.raises(ValueError) as error:
        locus.to_coordinate(Point(position=-5, offset=0))
    assert str(error.value) == 'Value must be non-negative.'

    with pytest.raises(IndexError) as error:
        locus.to_coordinate(Point(position=5, offset=0))
    assert str(error.value) == 'Position 5 exceeds locus length.'
    with pytest.raises(IndexError) as error:
        locus.to_coordinate(Point(position=0, offset=2))
    assert str(error.value) == 'Offset 2 should be at a locus end.'
    with pytest.raises(IndexError) as error:
        locus.to_coordinate(Point(position=4, offset=-2))
    assert str(error.value) == 'Offset -2 should be at a locus start.'
    with pytest.raises(ValueError) as error:
        locus.to_coordinate(Point(position=2, offset=1))
    assert str(error.value) == 'Position 2 is not at a locus boundary.'


def test_invalid_locus_inverted_point():
    """Reverse orientent Locus with invalid point."""
    locus = Locus((30, 35), True)

    with pytest.raises(ValueError) as error:
        locus.to_coordinate(Point(position=-5, offset=0))
    assert str(error.value) == 'Value must be non-negative.'
    with pytest.raises(IndexError) as error:
        locus.to_coordinate(Point(position=5, offset=0))
    assert str(error.value) == 'Position 5 exceeds locus length.'
    with pytest.raises(IndexError) as error:
        locus.to_coordinate(Point(position=0, offset=2))
    assert str(error.value) == 'Offset 2 should be at a locus end.'
    with pytest.raises(IndexError) as error:
        locus.to_coordinate(Point(position=4, offset=-2))
    assert str(error.value) == 'Offset -2 should be at a locus start.'
    with pytest.raises(ValueError) as error:
        locus.to_coordinate(Point(position=2, offset=1))
    assert str(error.value) == 'Position 2 is not at a locus boundary.'


def test_locus():
    """Forward orientent Locus."""
    locus = Locus((30, 35))

    invariant(locus.to_position, 29, locus.to_coordinate, Point(position=0, offset=-1))
    invariant(locus.to_position, 30, locus.to_coordinate, Point(position=0, offset=0))
    invariant(locus.to_position, 31, locus.to_coordinate, Point(position=1, offset=0))
    invariant(locus.to_position, 33, locus.to_coordinate, Point(position=3, offset=0))
    invariant(locus.to_position, 34, locus.to_coordinate, Point(position=4, offset=0))
    invariant(locus.to_position, 35, locus.to_coordinate, Point(position=4, offset=1))


def test_locus_inverted():
    """Reverse orientent Locus."""
    locus = Locus((30, 35), True)

    invariant(locus.to_position, 35, locus.to_coordinate, Point(position=0, offset=-1))
    invariant(locus.to_position, 34, locus.to_coordinate, Point(position=0, offset=0))
    invariant(locus.to_position, 33, locus.to_coordinate, Point(position=1, offset=0))
    invariant(locus.to_position, 31, locus.to_coordinate, Point(position=3, offset=0))
    invariant(locus.to_position, 30, locus.to_coordinate, Point(position=4, offset=0))
    invariant(locus.to_position, 29, locus.to_coordinate, Point(position=4, offset=1))
