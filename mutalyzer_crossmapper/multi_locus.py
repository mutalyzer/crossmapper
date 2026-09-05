from bisect import bisect_right
from itertools import accumulate
from dataclasses import dataclass

from .location import _nearest_location
from .locus import Locus, Coord, _check_locus, _check_positive_int
from .locus import Point as LocusPoint


@dataclass(slots=True)
class Point(LocusPoint):
    """Point dataclass"""

    region: str = ''

    def __post_init__(self) -> None:
        LocusPoint.__post_init__(self)
        if self.region not in ('', 'u', 'd'):
            raise ValueError(
                f'Region {self.region} is invalid, it must be a string from "", "u" or "d".'
            )


def _check_coordinate_within_length(value: int, length: int) -> None:
    """Check if a zero-based coordinate is within reference length."""
    if value >= length:
        raise ValueError(
            f'Coordinate {value} is not within the bounds of the reference length {length}.'
        )


def _check_location_end_within_length(value: int, length: int) -> None:
    """Check if a half-open interval end is within reference length."""
    if value > length:
        raise ValueError(
            f'Location end {value} is inconsistent with reference length {length}.'
        )


def _check_multi_locus(locations: list[tuple[int, int]], length: int | None = None) -> None:
    """Check if the locations list is valid."""
    if not locations or not isinstance(locations, list):
        raise ValueError('Locations must be a non-empty list of tuples.')

    for locus in locations:
        _check_locus(locus)

    for l1, l2 in zip(locations, locations[1:]):
        if l2[0] < l1[0]:
            raise ValueError(f'Locus {l1} and locus {l2} are not in ascending order.')
        if l2[0] < l1[1]:
            raise ValueError(f'Locus {l1} and locus {l2} are overlapping.')

    if length is not None:
        _check_positive_int(length)
        _check_location_end_within_length(locations[-1][1], length)


def _offsets(locations: list[tuple[int, int]], orientation: int) -> list[int]:
    """For each location, calculate the length of the preceding locations.

    :arg list locations: List of locations.
    :arg int orientation: Direction of locations.

    :returns list: List of cumulative location lengths.
    """
    return [0] + list(accumulate(map(lambda x: x[1] - x[0], locations[::orientation][:-1])))


class MultiLocus():
    """MultiLocus object."""

    def __init__(
        self,
        locations: list[tuple[int, int]],
        inverted: bool = False,
        length: int | None = None,
    ) -> None:
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        :arg int|None length: Length of the reference sequence, None if unknown.
        """
        _check_multi_locus(locations, length)
        if not isinstance(inverted, bool):
            raise ValueError(f'Value {inverted} is not a boolean.')
        self._locations = locations
        self._inverted = inverted
        self._length = length

        self._loci = [Locus(location, inverted) for location in locations]
        self._orientation = -1 if inverted else 1
        self._offsets = _offsets(locations, self._orientation)
        # one-based length of the MultiLocus
        self._end = sum(end - start for start, end in locations)

    def _validate_coord(self, coordinate: int) -> None:
        """Check if the coordinate is valid."""
        if self._length is not None:
            _check_coordinate_within_length(coordinate, self._length)

    def _validate_point(self, index: int, position: int, offset: int, region: str) -> None:
        """Validate if a multi locus Point is valid under HGVS rules.

        :arg int index: Index of the locus.
        :arg int position: Position.
        :arg int offset: Offset.
        :arg str region: Region.
        """
        if region == 'u':
            if position != 0:
                raise ValueError(f'Position {position} is not at upstream boundary.')
            if offset >= 0:
                raise ValueError(f'Offset {offset} at upstream region should be negative.')
            if self._inverted:
                if (
                    self._length is not None
                    and -offset
                    >= self._length - self._loci[self._direction(0)].boundary[1]
                ):
                    raise ValueError(f'Offset {offset} exceeds upstream region.')
            else:
                if -offset > self._loci[self._direction(0)].boundary[0]:
                    raise ValueError(f'Offset {offset} exceeds upstream region.')

        if region == 'd':
            if position != self._end - 1:
                raise ValueError(f'Position {position} is not at downstream boundary.')
            if offset <= 0:
                raise ValueError(f'Offset {offset} at downstream region should be positive.')
            if not self._inverted:
                if (
                    self._length is not None
                    and offset
                    >= self._length - self._loci[self._direction(-1)].boundary[1]
                ):
                    raise ValueError(f'Offset {offset} exceeds downstream region.')
            else:
                if (
                    offset
                    > self._loci[self._direction(len(self._locations) - 1)].boundary[0]
                ):
                    raise ValueError(f'Offset {offset} exceeds downstream region.')

        if region == '':
            if offset < 0:
                if index == 0 and position == 0:
                    raise ValueError(f'Offset {offset} at the first locus should be in the upstream region.')
                if self._inverted:
                    if (
                        self._direction(index) != len(self._loci) - 1
                        and -offset
                        >= self._loci[self._direction(index - 1)].boundary[0]
                        - self._loci[self._direction(index)].boundary[1]
                    ):
                        raise IndexError(f'Offset {offset} exceeds intron length.')
                else:
                    if (
                        self._direction(index) != 0
                        and -offset
                        >= self._loci[self._direction(index)].boundary[0]
                        - self._loci[self._direction(index - 1)].boundary[1]
                    ):
                        raise IndexError(f'Offset {offset} exceeds intron length.')
            if offset > 0:
                if index == len(self._loci) - 1 and position == self._end - 1:
                    raise ValueError(f'Offset {offset} at the last locus should be in the downstream region.')
                if self._inverted:
                    if (
                        self._direction(index) != 0
                        and offset
                        >= self._loci[self._direction(index)].boundary[0]
                        - self._loci[self._direction(index + 1)].boundary[1]
                    ):
                        raise IndexError(f'Offset {offset} exceeds intron length.')
                else:
                    if (
                        self._direction(index) != len(self._loci) - 1
                        and offset
                        >= self._loci[self._direction(index + 1)].boundary[0]
                        - self._loci[self._direction(index)].boundary[1]
                    ):
                        raise IndexError(f'Offset {offset} exceeds intron length.')

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
        """Convert a coordinate dataclass to a multi locus point dataclass.

        :arg Coord coord: Coordinate dataclass.

        :returns Point: Multi locus point dataclass.
        """
        self._validate_coord(coord.coordinate)
        index = _nearest_location(self._locations, coord.coordinate, self._inverted)
        outside = self._orientation * self._outside(coord.coordinate)
        region = 'u' if outside < 0 else 'd' if outside > 0 else ''
        point = self._loci[index].to_position(coord)

        return Point(
            position=point.position + self._offsets[self._direction(index)],
            offset=point.offset,
            region=region,
        )

    def to_coordinate(self, point: Point) -> Coord:
        """Convert a multi locus point dataclass to a coordinate dataclass.

        :arg Point point: Multi locus point dataclass.

        :returns Coord: Coordinate dataclass.
        """
        index = min(
            len(self._offsets), max(0, bisect_right(self._offsets, point.position) - 1)
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

        except ValueError as error:
            if 'Position' in str(error):
                raise ValueError(
                    str(error).replace(
                        str(point.position - self._offsets[index]), str(point.position)
                    )
                ) from error
            raise
        except IndexError as error:
            if 'Position' in str(error):
                raise IndexError(
                    f'Position {point.position} exceeds multi locus length.'
                ) from error
            raise
