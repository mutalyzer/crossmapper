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
    def __init__(self, locations: list[tuple[int, int]], inverted: bool=False) -> None:
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
        """Convert a coordinate to a position.

        :arg int coordinate: Coordinate.

        :returns dict: Position model 'position', 'offset' and 'region' keys.
        """
        index = nearest_location(self._locations, coordinate, self._inverted)
        outside = self._orientation * self.outside(coordinate)
        region = 'u' if outside < 0 else 'd' if outside > 0 else ''
        locus_pos_m = self._loci[index].to_position(coordinate)

        if outside:
            return {
                'position': abs(locus_pos_m['offset']) - 1,
                'offset': 0,
                'region': region
            }
        return {
            'position': locus_pos_m['position'] + self._offsets[self._direction(index)],
            'offset': locus_pos_m['offset'],
            'region': region
        }

    def to_coordinate(self, pos_m: dict[str, int | str]) -> int:
        """Convert a position model to a coordinate.

        :arg dict pos_m: Position model with 'position','offset' and 'region' keys.

        :returns int: Coordinate.
        """
        region = pos_m['region']

        if region == 'u':
            if self._inverted:
                return self._locations[-1][1] + abs(pos_m['position']) - pos_m['offset']
            return self._locations[0][0] - abs(pos_m['position']) + pos_m['offset'] - 1
        elif region == 'd':
            if self._inverted:
                return self._locations[0][0] - abs(pos_m['position']) - pos_m['offset'] - 1
            return self._locations[-1][1] + abs(pos_m['position']) + pos_m['offset']

        index = min(
            len(self._offsets),
            max(0, bisect_right(self._offsets, pos_m['position']) - 1)
        )
        locus_pos_m = {**pos_m, 'position': pos_m['position'] - self._offsets[index]}
        return self._loci[self._direction(index)].to_coordinate(locus_pos_m)
