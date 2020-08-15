from bisect import bisect_right
from itertools import accumulate

from .location import nearest_location
from .locus import Locus


def _offsets(locations, orientation):
    """For each location, calculate the length of the preceding locations.

    :arg list locations: List of locations.
    :arg int orientation: Direction of {locations}.

    :returns list: List of cumulative location lengths.
    """
    return [0] + list(accumulate(map(
        lambda x: x[1] - x[0], locations[::orientation][:-1])))


class MultiLocus(object):
    """MultiLocus object."""
    def __init__(self, locations, inverted=False):
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        self._locations = locations
        self._inverted = inverted

        self._loci = [Locus(location, inverted) for location in locations]
        self._orientation = -1 if inverted else 1
        self._offsets = _offsets(locations, self._orientation)

    def _direction(self, index):
        if self._inverted:
            return len(self._offsets) - index - 1
        return index

    def outside(self, coordinate):
        """Calculate the offset relative to this MultiLocus.

        :arg int coordinate: Coordinate.

        :returns int: Negative: upstream, 0: inside, positive: downstream.
        """
        if coordinate < self._loci[0].boundary[0]:
            return coordinate - self._loci[0].boundary[0]
        if coordinate > self._loci[-1].boundary[1]:
            return coordinate - self._loci[-1].boundary[1]
        return 0

    def to_position(self, coordinate):
        """Convert a coordinate to a position.

        :arg int coordinate: Coordinate.

        :returns tuple: Position.
        """
        index = nearest_location(self._locations, coordinate, self._inverted)
        outside = self._orientation * self.outside(coordinate)
        location = self._loci[index].to_position(coordinate)

        return (
            location[0] + self._offsets[self._direction(index)],
            location[1],
            outside)

    def to_coordinate(self, position):
        """Convert a position to a coordinate.

        :arg int position: Position.

        :returns int: Coordinate.
        """
        index = min(
            len(self._offsets),
            max(0, bisect_right(self._offsets, position[0]) - 1))

        return self._loci[self._direction(index)].to_coordinate(
            (position[0] - self._offsets[index], position[1]))
