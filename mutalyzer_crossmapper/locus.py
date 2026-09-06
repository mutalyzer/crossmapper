"""Conversions between coordinates and points for a single locus.

A locus is a half-open interval ``[start, end)`` defined on a sequence.
It can be addressed either by a coordinate, which is a zero-based index into
that sequence, or by a point, which is a zero-based position within the locus
plus an offset relative to its start or end.
"""
from dataclasses import dataclass


def _check_int(value: int) -> None:
    """Check if the input value type is integer."""
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError('Value must be an integer.')


def _check_non_negative_int(value: int) -> None:
    """Check if the value is a non-negative integer."""
    _check_int(value)
    if value < 0:
        raise ValueError('Value must be non-negative.')


def _check_positive_int(value: int) -> None:
    """Check if the value is a positive integer."""
    _check_int(value)
    if value < 1:
        raise ValueError(f'Value {value} is not positive.')


@dataclass(slots=True)
class Point:
    """A position within a locus, with an offset from its start or end.

    The offset is zero for positions inside the locus. A non-zero offset is
    only valid at a boundary: negative at the start, positive at the end.

    :arg int position: Zero-based position within the locus.
    :arg int offset: Offset relative to the locus start or end.
    """
    position: int
    offset: int = 0

    def __post_init__(self) -> None:
        _check_non_negative_int(self.position)
        _check_int(self.offset)


def _check_locus(location: tuple[int, int]) -> None:
    """Check if the locus location is valid."""
    if not isinstance(location, tuple) or len(location) != 2:
        raise ValueError('Locus must be a tuple of two values.')

    for value in location:
        _check_non_negative_int(value)

    if location[0] >= location[1]:
        raise ValueError(
            f'Locus start {location[0]} must be smaller than locus end {location[1]}.')


class Locus():
    """Convert coordinates to and from points within a single locus.

    For ``Locus((10, 15))``, the half-open interval contains coordinates 10 to
    14 and positions 0 to 4. Coordinates outside the locus are represented
    as boundary positions with offsets::

       coordinate        8    9   10   11   12   13   14   15   16
                         |    |    |    |    |    |    |    |    |
       point position    0    0    0    1    2    3    4    4    4
             offset     -2   -1    0    0    0    0    0    1    2

    >>> from mutalyzer_crossmapper.locus import Locus, Point
    >>> locus = Locus((10, 15))
    >>> locus.to_position(9)
    Point(position=0, offset=-1)
    >>> locus.to_position(15)
    Point(position=4, offset=1)
    >>> locus.to_coordinate(Point(position=4, offset=1))
    15

    For a locus on the reverse-complement strand, set ``inverted=True``::

       coordinate        8    9   10   11   12   13   14   15   16
                         |    |    |    |    |    |    |    |    |
       point position    4    4    4    3    2    1    0    0    0
             offset      2    1    0    0    0    0    0   -1   -2

    >>> inverted = Locus((10, 15), inverted=True)
    >>> inverted.to_position(14)
    Point(position=0, offset=0)
    >>> inverted.to_position(15)
    Point(position=0, offset=-1)
    """

    def __init__(self, location: tuple[int, int], inverted: bool = False) -> None:
        """Initialize a Locus object.

        :arg tuple location: Half-open interval, with the start strictly
            smaller than the end.
        :arg bool inverted: Orientation.

        :raises ValueError: If the location is not a tuple of two non-negative
            integers, or its start is not smaller than its end.
        """
        _check_locus(location)

        self._inverted = inverted
        self.boundary = location[0], location[1] - 1
        self._end = location[1] - location[0]  # one-based length of the locus

    def _validate_point(self, position: int, offset: int) -> None:
        """Check that a non-zero offset is at the matching locus boundary.

        Negative offsets belong at the start, positive ones at the end, and the
        position must lie within the locus.
        """
        if offset != 0 and position not in (0, self._end - 1):
            raise ValueError(f'Position {position} is not at a locus boundary.')
        if offset < 0 and position != 0:
            raise IndexError(f'Offset {offset} should be at a locus start.')
        if offset > 0 and position != self._end - 1:
            raise IndexError(f'Offset {offset} should be at a locus end.')
        if position > self._end - 1:
            raise IndexError(f'Position {position} exceeds locus length.')

    def to_position(self, coordinate: int) -> Point:
        """Convert a coordinate to a locus point dataclass.

        Coordinates outside the locus are converted to a boundary position
        with a non-zero offset.

        :arg int coordinate: Zero-based coordinate, non-negative.

        :returns Point: Locus point dataclass.

        :raises ValueError: If the coordinate is not a non-negative integer.
        """
        _check_non_negative_int(coordinate)

        if self._inverted:
            if coordinate > self.boundary[1]:
                return Point(position=0, offset=self.boundary[1] - coordinate)
            if coordinate < self.boundary[0]:
                return Point(position=self._end - 1, offset=self.boundary[0] - coordinate)
            return Point(position=self.boundary[1] - coordinate, offset=0)

        if coordinate < self.boundary[0]:
            return Point(position=0, offset=coordinate - self.boundary[0])
        if coordinate > self.boundary[1]:
            return Point(position=self._end - 1, offset=coordinate - self.boundary[1])
        return Point(position=coordinate - self.boundary[0], offset=0)

    def to_coordinate(self, point: Point) -> int:
        """Convert a locus point dataclass to a coordinate.

        :arg Point point: Locus point dataclass, its position relative to this
            locus.

        :returns int: Coordinate.

        :raises ValueError: If a non-zero offset is not at a locus boundary, or
            the point converts to a negative coordinate.
        :raises IndexError: If the position exceeds the locus length, or the
            offset is at the wrong boundary.
        """
        self._validate_point(point.position, point.offset)

        if self._inverted:
            coordinate = self.boundary[1] - point.position - point.offset
        else:
            coordinate = self.boundary[0] + point.position + point.offset
        if coordinate < 0:
            raise ValueError(
                f'Position {point.position} with offset {point.offset} converts to '
                f'negative coordinate {coordinate}.'
            )

        return coordinate
