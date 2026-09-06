from dataclasses import dataclass


@dataclass(slots=True)
class Point:
    """Locus point dataclass."""
    position: int
    offset: int = 0

    def __post_init__(self) -> None:
        _check_non_negative_int(self.position)
        _check_int(self.offset)


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


def _check_locus(locus: tuple[int, int]) -> None:
    """Check if the locus location is valid."""
    if not isinstance(locus, tuple) or len(locus) != 2:
        raise ValueError('Locus must be a tuple of two values.')

    for value in locus:
        _check_non_negative_int(value)

    if locus[0] >= locus[1]:
        raise ValueError(f'Locus start {locus[0]} must be smaller than locus end {locus[1]}.')


class Locus():
    """Locus object."""

    def __init__(self, location: tuple[int, int], inverted: bool = False) -> None:
        """Initialize a Locus object.

        :arg tuple location: Locus location.
        :arg bool inverted: Orientation.
        """
        _check_locus(location)

        self._inverted = inverted
        self.boundary = location[0], location[1] - 1
        self._end = location[1] - location[0]  # one-based length of the locus

    def _validate_point(self, position: int, offset: int) -> None:
        """Validate a locus point dataclass according to HGVS rules."""
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

        :arg int coordinate: Coordinate.

        :returns Point: Locus point dataclass.
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

        :arg Point point: Locus point dataclass.

        :returns int: Coordinate.
        """
        self._validate_point(point.position, point.offset)

        if self._inverted:
            return self.boundary[1] - point.position - point.offset
        return self.boundary[0] + point.position + point.offset
