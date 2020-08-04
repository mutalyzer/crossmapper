"""Crossmapper position conversion library.

Definitions:

- Coordinates are zero based, non-negative integers.
- Locations are zero based right-open non-negative integer intervals,
  consistent with Python's range() and sequence slicing functions.
- Loci and exons are locations.
- An exon list is a list of locations that, when flattened, is an increasing
  sequence.
- A position is a 2-tuple of which the first element is a one based non-zero
  integer relative to an element in a location and the second element is an
  integer offset relative to the first element.
"""
from bisect import bisect_left

from .location import cut_locations, nearest_location


def _offsets(locations, inverted=False):
    """For each location, calculate the length of the preceding locations.

    :arg list locations: List of locations.
    :arg bool inverted: Direction of {locations}.

    :returns list: List of cumulative location lengths.
    """
    lengths = []

    length = 0
    direction = -1 if inverted else 1

    for location in locations[::direction]:
        lengths.append(length)
        length += location[1] - location[0]

    return lengths


class Locus(object):
    """Locus object."""
    def __init__(self, location, inverted=False):
        """
        :arg tuple location: Locus location.
        :arg bool inverted: Orientation.
        """
        self._location = location
        self._inverted = inverted

    def _to_position(self, coordinate):
        """Convert a coordinate to a proper position.

        :arg int coordinate: Coordinate.

        :returns tuple: Position.
        """
        if self._inverted:
            if coordinate >= self._location[1]:
                return 1, self._location[1] - coordinate - 1
            if coordinate < self._location[0]:
                return (
                    self._location[1] - self._location[0],
                    self._location[0] - coordinate)
            return self._location[1] - coordinate, 0

        if coordinate < self._location[0]:
            return 1, coordinate - self._location[0]
        if coordinate >= self._location[1]:
            return (
                self._location[1] - self._location[0],
                coordinate - self._location[1] + 1)
        return coordinate - self._location[0] + 1, 0

    def to_position(self, coordinate, degenerate=False):
        """Convert a coordinate to a position.

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns tuple: Position.
        """
        position = self._to_position(coordinate)

        if degenerate:
            if position[0] == 1:
                return position[1], 0
            return position[0] + position[1], 0

        return position

    def to_coordinate(self, position):
        """Convert a position to a coordinate.

        :arg int position: Position.

        :returns int: Coordinate.
        """
        if self._inverted:
            if position[0] < 1:  # Degenerate position.
                return self._location[1] - position[0] - position[1] - 1
            return self._location[1] - position[0] - position[1]

        if position[0] < 1:      # Degenerate position.
            return self._location[0] + position[0] + position[1]
        return self._location[0] + position[0] + position[1] - 1


class MultiLocus(object):
    """MultiLocus object."""
    def __init__(self, locations, inverted=False, negated=False):
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        :arg bool negated: Change the sign of all positions.
        """
        self._locations = locations
        self._inverted = inverted
        self._negated = negated

        self._loci = [Locus(location, inverted) for location in locations]
        self._offsets = _offsets(locations, inverted)

        self.orientation = 1
        if self._inverted != self._negated:
            self.orientation = -1

    def _sign(self, position):
        if self._negated:
            return -position[0], -position[1], position[2]
        return position

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

        return self._sign((
            location[0] + self._offsets[self._direction(index)],
            location[1],
            outside))

    def to_coordinate(self, position):
        """Convert a position to a coordinate.

        :arg int position: Position.

        :returns int: Coordinate.
        """
        position_ = self._sign(position)

        index = min(
            len(self._offsets),
            max(0, bisect_left(self._offsets, position_[0]) - 1))

        return self._loci[self._direction(index)].to_coordinate(
            (position_[0] - self._offsets[index], position_[1]))


class Crossmap(object):
    """Crossmap object."""
    _noncoding_error = 'no locations provided'
    _coding_error = 'no cds provided'

    def __init__(self, locations=None, cds=None, inverted=False):
        """
        :arg list locations: List of locus locations.
        :arg tuple cds: Locus location.
        :arg bool inverted: Orientation.
        """
        self._inverted = inverted

        self._noncoding = None
        if locations:
            self._noncoding = MultiLocus(locations, inverted)

        self._coding = None
        if cds:
            head, tail = cut_locations(locations, cds[0])
            regions = (head, *cut_locations(tail, cds[1]))

            self._coding = (
                MultiLocus(regions[0], True, not inverted),
                MultiLocus(regions[1], inverted),
                MultiLocus(regions[2], False, inverted))

            self._regions = [(x[0][0], x[-1][1]) if x else () for x in regions]
            self._cds_len = sum(map(lambda x: x[1] - x[0], regions[1]))

    def _check(self, condition, error):
        if not condition:
            raise ValueError(error)

    def _nearest_region(self, coordinate):
        """
        :returns int: -1: left, 0: middle, 1: right.
        """
        outside = self._coding[1].outside(coordinate)

        if outside < 0 and self._regions[0]:
            return nearest_location(
                self._regions[:2], coordinate, self._inverted) - 1
        if outside > 0 and self._regions[2]:
            return nearest_location(
                self._regions[1:], coordinate, self._inverted)

        return 0

    def coordinate_to_genomic(self, coordinate):
        """Convert a coordinate to a genomic position (g./m./o.).

        :arg int coordinate: Coordinate.

        :returns int: Genomic position.
        """
        return coordinate + 1

    def genomic_to_coordinate(self, position):
        """Convert a genomic position (g./m./o.) to a coordinate.

        :arg int position: Genomic position.

        :returns int: Coordinate.
        """
        return position - 1

    def coordinate_to_noncoding(self, coordinate):
        """Convert a coordinate to a noncoding position (n./r.).

        :arg int coordinate: Coordinate.

        :returns tuple: Noncoding position.
        """
        self._check(self._noncoding, self._noncoding_error)

        return self._noncoding.to_position(coordinate)

    def noncoding_to_coordinate(self, position):
        """Convert a noncoding position (n./r.) to a coordinate.

        :arg tuple position: Noncoding position.

        :returns int: Coordinate.
        """
        self._check(self._noncoding, self._noncoding_error)

        return self._noncoding.to_coordinate(position)

    def coordinate_to_coding(self, coordinate, degenerate=False):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns tuple: Coding position (c./r.).
        """
        self._check(self._coding, self._coding_error)

        region = self._nearest_region(coordinate)
        position = self._coding[region + 1].to_position(coordinate)
        #selected_region = self._coding[1].orientation * region + 1#_region(region, self._inverted)

        #r = self._noncoding.to_position(coordinate)[2]
        #if degenerate and r != 1:
        #    if selected_region == 1:
        #        if r == 0 and position[0] == 1 and position[1] < 0:
        #            return (position[1], 0, 0)
        #        if r == 2 and position[0] == self._cds_len and position[1] > 0:
        #            return (position[1], 0, 2)
        #    return (position[0] + position[1], 0, selected_region)

        outside = self._coding[1].orientation * region

        # Remove "outside" coordinates between regions.
        if position[2] * outside < 0:
            return (position[0], position[1], 0, outside)
        if not outside:
            if position[2] < 0 and self._regions[0]:
                return (position[0], position[1], 0, outside)
            if position[2] > 0 and self._regions[2]:
                return (position[0], position[1], 0, outside)

        return (*position, outside)

    def coding_to_coordinate(self, position):
        """Convert a coding position (c./r.) to a coordinate.

        :arg tuple position: Coding position (c./r.).

        :returns int: Coordinate.
        """
        self._check(self._coding, self._coding_error)

        region = self._coding[1].orientation * position[3] + 1
        #region = _region(position[2], self._inverted)
        #if not self._regions[region]:  # Degenerate position.
        #    if position[2] == 0:
        #        return self._coding[1].to_coordinate(position[:2])
        #    if position[2] == 2:
        #        return self._coding[1].to_coordinate(
        #            (self._cds_len + position[0], position[1], 1))

        return self._coding[region].to_coordinate((*position[:2], 0))

    def coordinate_to_protein(self, coordinate):
        """Convert a coordinate to a protein position (p.).

        :arg int coordinate: Coordinate.

        :returns tuple: Protein position (p.).
        """
        self._check(self._coding, self._coding_error)

        position = self.coordinate_to_coding(coordinate)

        if not position[2]:
            return position[0] // 3, position[0] % 3 + 1, position[1], 0
        return (
            (position[0] + 2) // 3, (position[0] + 2) % 3 + 1,
            position[1], position[2])

    def protein_to_coordinate(self, position):
        """Convert a protein position (p.) to a coordinate.

        :arg tuple position: Protein position (p.).

        :returns int: Coordinate.
        """
        self._check(self._coding, self._coding_error)

        if not position[3]:
            return self.coding_to_coordinate(
                (position[0] * 3 + position[1] - 1, position[2], position[3]))

        return self.coding_to_coordinate(
            (position[0] * 3 - 3 + position[1], position[2], position[3]))
