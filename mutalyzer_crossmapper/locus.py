from dataclasses import dataclass
from .checker import _check_locus, _check_non_negative, _check_int


@dataclass(slots=True)
class Point:
    """Point dataclass"""
    position: int
    offset: int = 0

    def __post_init__(self) -> None:
        _check_non_negative(self.position)
        _check_int(self.offset)


@dataclass(slots=True)
class Coord:
    """Coordinate dataclass"""
    coordinate: int

    def __post_init__(self) -> None:
        _check_non_negative(self.coordinate)


class Locus(object):
    """Locus object."""

    def __init__(self, location: tuple[int, int], inverted: bool = False) -> None:
        """
        :arg tuple location: Locus location.
        :arg bool inverted: Orientation.
        """
        _check_locus(location)

        self._inverted = inverted
        self.boundary = location[0], location[1] - 1  #: 0 based open on coordinate
        self._end = self.boundary[1] - self.boundary[0]

    def _validate_point(self, position, offset) -> None:
        """Validate a point model under HGVS rules.

        :arg int position: Position.
        :arg int offset: Offset.
        """
        if offset != 0 and position not in (0, self._end):
            raise ValueError(f"Position {position} is not at locus boundary.")
        if offset < 0 and position != 0:
            raise IndexError(f"Offset {offset} should be at a locus start.")
        if offset > 0 and position != self._end:
            raise IndexError(f"Offset {offset} should be at a locus end.")
        if position > self._end:
            raise IndexError(f"Position {position} exceeds locus length {self._end + 1}")

    def to_position(self, coord: Coord) -> Point:
        """Convert a coordinate to a proper point model.

        :arg Coord coord: Coordinate module.

        :returns Point: Position point model.
        """
        if self._inverted:
            if coord.coordinate > self.boundary[1]:
                return Point(position=0, offset=self.boundary[1] - coord.coordinate)
            if coord.coordinate < self.boundary[0]:
                return Point(position=self._end, offset=self.boundary[0] - coord.coordinate)
            return Point(position=self.boundary[1] - coord.coordinate, offset=0)

        if coord.coordinate < self.boundary[0]:
            return Point(position=0, offset=coord.coordinate - self.boundary[0])
        if coord.coordinate > self.boundary[1]:
            return Point(position=self._end, offset=coord.coordinate - self.boundary[1])
        return Point(position=coord.coordinate - self.boundary[0], offset=0)

    def to_coordinate(self, point: Point) -> Coord:
        """Convert a point model to a coordinate model.

        :arg Point point: Point model.

        :returns Coord: Coordinate model.
        """
        self._validate_point(point.position, point.offset)

        if self._inverted:
            return Coord(coordinate=self.boundary[1] - point.position - point.offset)
        return Coord(coordinate=self.boundary[0] + point.position + point.offset)
