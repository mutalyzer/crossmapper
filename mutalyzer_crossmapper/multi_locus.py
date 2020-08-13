from bisect import bisect_left

from .location import nearest_location
from .locus import Locus


def _offsets(locations, inverted=False):
    """For each location, calculate the length of the preceding locations.

    :arg list locations: List of locations.
    :arg int orientation: Direction of {locations}.

    :returns list: List of cumulative location lengths.
    """
    lengths = []

    length = 0
    direction = -1 if inverted else 1

    for location in locations[::direction]:
        lengths.append(length)
        length += location[1] - location[0]

    return lengths


class MultiLocus(object):
    """MultiLocus object."""
    #def __init__(self, locations, inverted=False, negated=False):
    def __init__(self, locations, inverted=False):
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        #:arg bool negated: Change the sign of all positions.
        """
        self._locations = locations
        self._inverted = inverted
        #self._negated = negated

        self._loci = [Locus(location, inverted) for location in locations]
        self._offsets = _offsets(locations, inverted)

        self.orientation = 1
        #if self._inverted != self._negated:
        if self._inverted:
            self.orientation = -1

    #def _sign(self, position):
    #    if self._negated:
    #        return -position[0], -position[1], position[2]
    #    return position

    def _direction(self, index):
        if self._inverted:
            return len(self._offsets) - index - 1
        return index

    def outside(self, coordinate):
        """Calculate the offset relative to this MultiLocus.

        :arg int coordinate: Coordinate.

        :returns int: Negative: upstream, 0: inside, positive: downstream.
        """
        if coordinate < self._locations[0][0]:
            return coordinate - self._locations[0][0]
        if coordinate >= self._locations[-1][1]:
            return coordinate - self._locations[-1][1] + 1
        return 0

    def to_position(self, coordinate, degenerate=False):
        """Convert a coordinate to a position.

        :arg int coordinate: Coordinate.

        :returns tuple: Position.
        """
        index = nearest_location(self._locations, coordinate, self._inverted)
        outside = self.orientation * self.outside(coordinate)
        location = self._loci[index].to_position(
            coordinate, outside and degenerate)

        return (
            location[0] + self._offsets[self._direction(index)],
            location[1],
            outside)
        #return self._sign((
        #    location[0] + self._offsets[self._direction(index)],
        #    location[1],
        #    outside))

    def to_coordinate(self, position):
        """Convert a position to a coordinate.

        :arg int position: Position.

        :returns int: Coordinate.
        """
        #position_ = self._sign(position)

        #index = min(
        #    len(self._offsets),
        #    max(0, bisect_left(self._offsets, position_[0]) - 1))

        #return self._loci[self._direction(index)].to_coordinate(
        #    (position_[0] - self._offsets[index], position_[1]))
        index = min(
            len(self._offsets),
            max(0, bisect_left(self._offsets, position[0]) - 1))

        return self._loci[self._direction(index)].to_coordinate(
            (position[0] - self._offsets[index], position[1]))
