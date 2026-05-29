from bisect import bisect_right
from itertools import accumulate

from .location import nearest_location
from .locus import Locus


def _offsets(locations: list[tuple[int, int]], orientation: int) -> list[int]:
    """For each location, calculate the length of the preceding locations.

    :arg list locations: List of locations.
    :arg int orientation: Direction of {locations}.

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
        self._locations = locations
        self._inverted = inverted

        self._loci = [Locus(location, inverted) for location in locations]
        self._orientation = -1 if inverted else 1
        self._offsets = _offsets(locations, self._orientation)

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

    def to_position(self, coordinate: int) -> dict[str, int | str]:
        """Convert a coordinate to a point model.

        :arg int coordinate: Coordinate.

        :returns dict: Point model 'position', 'offset' and 'region' keys.
        """
        index = nearest_location(self._locations, coordinate, self._inverted)
        outside = self._orientation * self.outside(coordinate)
        region = 'u' if outside < 0 else 'd' if outside > 0 else ''
        point = self._loci[index].to_position(coordinate)
        return {
            'position': point['position'] + self._offsets[self._direction(index)],
            'offset': point['offset'],
            'region': region,
        }

    def to_coordinate(self, point: dict[str, int | str]) -> int:
        """Convert a point model to a coordinate.

        :arg dict point: Point model with 'position','offset' and 'region' keys.

        :returns int: Coordinate.
        """
        region = point['region']

        if region == 'u':
            if self._inverted:
                return self._locations[-1][1] - point['offset'] - 1
            return self._locations[0][0] + point['offset']
        if region == 'd':
            if self._inverted:
                return self._locations[0][0] - point['offset']
            return self._locations[-1][1] + point['offset'] - 1

        index = min(
            len(self._offsets),
            max(0, bisect_right(self._offsets, point['position']) - 1)
        )
        return self._loci[self._direction(index)].to_coordinate(
            {
                'position': point['position'] - self._offsets[index],
                'offset': point['offset'],
                'region': point['region'],
            }
        )
