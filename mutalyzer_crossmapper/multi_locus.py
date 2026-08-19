from bisect import bisect_right
from itertools import accumulate
from dataclasses import dataclass


from .location import nearest_location
from .locus import Locus, Coord, Point as LocusPoint, _check_locus


@dataclass(slots=True)
class Point(LocusPoint):
    """Point dataclass"""
    region: str = ''

    def __post_init__(self) -> None:
        LocusPoint.__post_init__(self)
        if self.region not in ('', 'u', 'd'):
            raise ValueError(f"Region {self.region} is not valid. Must be '', 'u', or 'd'.")


def _check_in_range(value: int, length: int | None = None) -> None:
    """Check if the value no larger than length."""
    print(value, length)
    if length is not None and value >= length:
        raise ValueError(f"Value {value} must be within the bounds of the reference sequence {length}.")


def _check_multi_locus(locations: list[tuple[int, int]], length: int | None = None) -> None:
    print(locations, length)
    """Check if the locations list is valid."""
    for locus in locations:
        _check_locus(locus)

    for l1, l2 in zip(locations, locations[1:]):
        if l2[0] < l1[1]:
            raise ValueError(f"Locus {l2} and locus {l1} are overlapping.")

    _check_in_range(locations[-1][1], length)


def _offsets(locations: list[tuple[int, int]], orientation: int) -> list[int]:
    """For each location, calculate the length of the preceding locations.

    :arg list locations: List of locations.
    :arg int orientation: Direction of locations.

    :returns list: List of cumulative location lengths.
    """
    return [0] + list(accumulate(map(
        lambda x: x[1] - x[0], locations[::orientation][:-1])))


class MultiLocus(object):
    """MultiLocus object."""
    def __init__(self, locations: list[tuple[int, int]], length: int |None = None, inverted: bool = False) -> None:
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        _check_multi_locus(locations, length)
        self._locations = locations
        self._inverted = inverted
        self._length = length

        self._loci = [Locus(location, inverted) for location in locations]
        self._orientation = -1 if inverted else 1
        self._offsets = _offsets(locations, self._orientation)
        self._end = sum(end - start for start, end in locations)    # one-based length of the MultiLocus

    def _validate_coord(self, coord) -> None:
        """Check if the coordinate is valid."""
        if self._length is not None:
            _check_in_range(coord, self._length)

    def _validate_point(self, index: int, position: int, offset: int, region: str) -> None:
        """Check if a point is valid under HGVS rules.

        :arg int index: Index of the locus.
        :arg int position: Position.
        :arg int offset: Offset.
        :arg str region: Region.
        """
        print("offsets", self._offsets)
        # Upstream region validation, position is constant value and offset should not be positive
        if region == 'u':
            if position != self._offsets[0]:
                raise ValueError(f"Position {position} is not at the upstream boundary.")
            if offset > 0:
                raise ValueError(f"Offset {offset} at upstream boundary should not be positive.")
            if self._inverted:
                if self._length is not None and abs(offset) >= self._length - self._loci[self._direction(0)].boundary[1]:
                    raise ValueError(f"Offset {offset} exceeds upstream region.")
            else:
                if abs(offset) > self._loci[self._direction(0)].boundary[0]:
                    raise ValueError(f"Offset {offset} exceeds upstream boundary.")

        # Downstream region validation, position is constant value and offset should not be negative
        if region == 'd':
            if position != self._end-1:
                raise ValueError(f"Position {position} is not at the downstream boundary {self._end-1}.")
            if offset < 0:
                raise ValueError(f"Offset {offset} at downstream boundary should not be negative.")
            if not self._inverted:
                if self._length is not None and abs(offset) >= self._length - self._loci[self._direction(-1)].boundary[1]:
                    raise ValueError(f"Offset {offset} exceeds downstream region.")
            else:
                if abs(offset) > self._loci[self._direction(0)].boundary[0]:
                    raise ValueError(f"Offset {offset} exceeds downstream boundary.")

        # '' region validation, position should be within the MultiLocus and offset should not exceed intron length
        if region == '':
            if position > self._end-1:
                raise IndexError(f"Position {position} exceeds MultiLocus length {self._end}")
            if offset < 0 and abs(offset) > abs(self._loci[self._direction(index)].boundary[0] - self._loci[self._direction(index-1)].boundary[1]):
                raise IndexError(f"Offset {offset} exceeds intron length.")
            if offset > 0 and abs(offset) > abs(self._loci[self._direction(index+1)].boundary[0] - self._loci[self._direction(index)].boundary[1]):
                raise IndexError(f"Offset {offset} exceeds intron length.")


    def _direction(self, index: int) -> int:
        if self._inverted:
            return len(self._offsets) - index - 1
        return index

    def _outside(self, coordinate: int) -> int:
        """Calculate the offset relative to this MultiLocus.

        :arg int coordinate: Coordinate.

        :returns int: Negative: upstream, 0: inside, positive: downstream.
        """
        if coordinate < self._loci[0].boundary[0]:
            return coordinate - self._loci[0].boundary[0]
        if coordinate > self._loci[-1].boundary[1]:
            return coordinate - self._loci[-1].boundary[1]
        return 0

    def to_position(self, coord: Coord) -> Point:
        """Convert a coordinate to a point model.

        :arg Coord coord: Coordinate model.

        :returns Point: Point model.
        """
        self._validate_coord(coord.coordinate)
        index = nearest_location(self._locations, coord.coordinate, self._inverted)
        outside = self._orientation * self._outside(coord.coordinate)
        region = 'u' if outside < 0 else 'd' if outside > 0 else ''
        point = self._loci[index].to_position(coord)

        return Point(
            position=point.position + self._offsets[self._direction(index)],
            offset=point.offset,
            region=region,
        )

    def to_coordinate(self, point: Point) -> int:
        """Convert a point model to a coordinate.

        :arg Point point: Point model.

        :returns int: Coordinate.
        """
        index = min(
            len(self._offsets),
            max(0, bisect_right(self._offsets, point.position) - 1)
        )
        self._validate_point(index, point.position, point.offset, point.region)

        if point.region == 'u':
            if self._inverted:
                return Coord(self._locations[-1][1] - point.offset - 1)
            return Coord(self._locations[0][0] + point.offset)
        if point.region == 'd':
            if self._inverted:
                return Coord(self._locations[0][0] - point.offset)
            return Coord(self._locations[-1][1] + point.offset - 1)

        try:
            return self._loci[self._direction(index)].to_coordinate(
                Point(
                    position=point.position - self._offsets[index],
                    offset=point.offset,
                )
            )

        except ValueError as e:
            if "Position" in str(e):
                raise ValueError(f"Position {point.position} is not at a locus boundary.") from e
            raise e
        except IndexError as e:
            if "Position" in str(e):
                raise IndexError(
                    f"Position {point.position} exceeds multi_locus length {self._end}"
                ) from e
            raise e

