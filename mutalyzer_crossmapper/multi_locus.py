from bisect import bisect_right
from itertools import accumulate
from operator import index

from .location import nearest_location
from .locus import Locus, Point


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
    def __init__(self, locations: list[tuple[int, int]], inverted: bool = False) -> None:
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        # Check if the locations are non-overlapping , e.g., [(1, 5), (4, 10)] should be invalid;
        # and sorted e.g., [(1, 5), (10, 15), (5, 10)] should be invalid.
        # Look for circular chromosome sequence
        self._locations = locations
        self._inverted = inverted

        self._loci = [Locus(location, inverted) for location in locations]
        self._orientation = -1 if inverted else 1
        self._offsets = _offsets(locations, self._orientation)

    # Consider add the length of the sequnce to the MultiLocus object,
    # so that we can check if a coordinate is outside the sequence length.
    def _validate_coordinate(self, coordinate: int) -> None:
        """Check if a coordinate is within the MultiLocus.

        :arg int coordinate: Coordinate.

        :raises ValueError: If coordinate is outside the MultiLocus.
        """
        if coordinate < 0:
            raise IndexError("Coordinate is outside sequence length.")

    def _validate_point(self, index: int, point: Point) -> None:
        """Check if a point is valid.

        :arg int index: Index of the locus.
        :arg Point point: Point model.

        :raises ValueError: If point is outside the MultiLocus.
        """
        if index == 0 and abs(point.offset) > self._loci[0].boundary[0]:
            raise IndexError(f"Offset {point.offset} is outside the intron length {self._loci[0].boundary[0]}.")
        if index > 0 and abs(point.offset) > self._loci[index].boundary[0] - self._loci[index - 1].boundary[1]:
            raise IndexError(f"Offset {point.offset} is outside the intron length {self._loci[index].boundary[0] - self._loci[index - 1].boundary[1]}.")

        if point.offset < 0:
            if point.position not in self._loci[index].boundary:
                raise ValueError(f"Position {point.position} is not at an exon boundary.")
            if self._loci[self._direction(index)].boundary[0] != point.position:
                raise IndexError(f"Offset {point.offset} should be '-' when position is at exon start.")

        if point.offset > 0:
            if point.position not in self._loci[index].boundary:
                raise ValueError(f"Position {point.position} is not at an exon boundary.")
            if self._loci[self._direction(index)].boundary[1] != point.position:
                raise IndexError(f"Offset {point.offset} should be '+' when position is at exon end.")


    def _direction(self, index: int) -> int:
        if self._inverted:
            return len(self._offsets) - index - 1
        return index

    def outside(self, coordinate: int) -> int:
        """Calculate the offset relative to this MultiLocus.

        :arg int coordinate: Coordinate.

        :returns int: Negative: upstream, 0: inside, positive: downstream.
        """
        if coordinate < self._loci[0].boundary[0]:
            return coordinate - self._loci[0].boundary[0]
        if coordinate > self._loci[-1].boundary[1]:
            return coordinate - self._loci[-1].boundary[1]
        return 0

    def to_position(self, coordinate: int) -> Point:
        """Convert a coordinate to a point model.

        :arg int coordinate: Coordinate.

        :returns Point: Point model .
        """
        self._validate_coordinate(coordinate)
        index = nearest_location(self._locations, coordinate, self._inverted)
        outside = self._orientation * self.outside(coordinate)
        region = 'u' if outside < 0 else 'd' if outside > 0 else ''
        point = self._loci[index].to_position(coordinate)
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


        if point.region == 'u':
            if self._inverted:
                return self._locations[-1][1] - point.offset - 1
            self._validate_point(0, point)
            return self._locations[0][0] + point.offset
        if point.region == 'd':
            if self._inverted:
                return self._locations[0][0] - point.offset
            return self._locations[-1][1] + point.offset - 1

        index = min(
            len(self._offsets),
            max(0, bisect_right(self._offsets, point.position) - 1)
        )
        self._validate_point(index, point)
        return self._loci[self._direction(index)].to_coordinate(
            Point(
                position=point.position - self._offsets[index],
                offset=point.offset,
                region=point.region,
            )
        )
