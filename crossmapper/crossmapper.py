"""Crossmapper position conversion library.

Definitions:

- Coordinates are zero based, non-negative integers.
- Ranges are zero based right-open non-negative integer intervals, consistent
  with Python's range() and sequence slicing functions.
- Loci and exons are ranges.
- An exon list is a list of ranges that, when flattened, is an increasing
  sequence.
- A position is a 2-tuple of which the first element is a one based non-zero
  integer relative to an element in a range and the second element is an
  integer offset relative to the first element.
"""
from bisect import bisect_left


def _nearest_boundary(coordinate, boundaries):
    """Given a coordinate, find the index of the nearest boundary.

    On a draw, the left boundary is chosen.

    :arg int coordinate: Coordinate.
    :arg list boundaries: List of boundaries.

    :returns int: Index of nearest boundary.
    """
    insertion_point = bisect_left(boundaries, coordinate)

    if (
            abs(coordinate - boundaries[insertion_point - 1]) <=
            abs(boundaries[insertion_point % len(boundaries)] - coordinate)):
        return insertion_point - 1
    return insertion_point


def _nearest_location(coordinate, boundaries):
    return _nearest_boundary(coordinate, boundaries) // 2


def _cut(coordinate, locations):
    """Divide a list of locations, cutting one of the locations in two.

    :arg int coordinate: Coordinate.
    :arg list locations: List of locations.

    :returns tuple: {locations} before coordinate, {locations} after coordinate.
    """
    boundaries = sum( # TODO: Fix this.
            [[location[0], location[1] - 1] for location in locations], [])

    location = _nearest_location(coordinate, boundaries)

    return (
        locations[:location] + [(locations[location][0], coordinate)],
        [(coordinate, locations[location][1])] + locations[location + 1:])


def _offsets(locations, inverted=False):
    """For each location, calculate the length of the preceding locations.

    :arg list locations: List of locations.
    :arg bool inverted: Direction of {locations}.

    :returns list: List of cumulative location lengths.
    """
    direction = 1
    if inverted:
        direction = -1

    length = 0
    lengths = []
    for location in locations[::direction]:
        lengths.append(length)
        length += location[1] - location[0]

    return lengths


class Locus(object):
    def __init__(self, location, inverted=False):
        """Construct a Locus object.

        :arg tuple location: Locus coordinates.
        :arg bool inverted: Orientation.
        """
        self._location = location
        self._inverted = inverted

    def to_position(self, coordinate):
        """Convert a coordinate to a position.

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

    def to_coordinate(self, position):
        """Convert a position to a coordinate.

        :arg int position: Position.

        :returns int: Coordinate.
        """
        if self._inverted:
            return self._location[1] - position[0] - position[1]
        return self._location[0] + position[0] + position[1] - 1


class MultiLocus(object):
    def __init__(self, locations, inverted=False, negated=False):
        """Construct a MultiLocus object.
        :arg list locations: List of locus coordinates.
        :arg bool inverted: Orientation.
        :arg bool negated: Change the sign of all positions.
        """
        self._inverted = inverted
        self._negated = negated

        self._loci = [Locus(location, inverted) for location in locations]
        self._boundaries = sum(
            [[location[0], location[1] - 1] for location in locations], [])
        self._offsets = _offsets(locations, inverted)

    def _sign(self, position):
        if self._negated:
            return -position[0], -position[1]
        return position

    def _direction(self, index):
        if self._inverted:
            return len(self._offsets) - index - 1
        return index

    def to_position(self, coordinate):
        """Convert a coordinate to a position.

        :arg int coordinate: Coordinate.

        :returns tuple: Position.
        """
        index = _nearest_location(coordinate, self._boundaries)
        location = self._loci[index].to_position(coordinate)

        return self._sign(
            (location[0] + self._offsets[self._direction(index)], location[1]))

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
    def __init__(self, locations=None, cds=None, inverted=False):
        """
        """
        self._cds = cds

        self._noncoding = None
        if locations:
            self._noncoding = MultiLocus(locations, inverted)

        self._coding = None
        if cds:
            utr5, tail = _cut(cds[0], locations)
            coding, utr3 = _cut(cds[1], tail)

            self._coding = (
                MultiLocus(utr5, not inverted, True),
                MultiLocus(coding, inverted),
                MultiLocus(utr3, inverted))

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
        return self._noncoding.to_position(coordinate)

    def noncoding_to_coordinate(self, position):
        """Convert a noncoding position (n./r.) to a coordinate.

        :arg tuple position: Noncoding position.

        :returns int: Coordinate.
        """
        return self._noncoding.to_coordinate(position)

    def coordinate_to_coding(self, coordinate):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.

        :returns tuple: Coding position (c./r.).
        """
        if coordinate < self._cds[0]:
            return (*self._coding[0].to_position(coordinate), 0)
        if coordinate >= self._cds[1]:
            return (*self._coding[2].to_position(coordinate), 2)
        return (*self._coding[1].to_position(coordinate), 1)

    def coding_to_coordinate(self, position):
        """Convert a coding position (c./r.) to a coordinate.

        :arg tuple position: Coding position (c./r.).

        :returns int: Coordinate.
        """
        return self._coding[position[2]].to_coordinate(position[:2])

#    def coordinate_to_protein(self, coordinate):
#        """Convert a coordinate to a protein position (p.).
#
#        Note that the converse of this function does not exist.
#
#        :arg int coordinate: Coordinate.
#
#        :returns tuple: Protein position (p.).
#        """
#        if not self._cds:
#            raise ValueError(
#                "conversion to protein position using a non coding transcript")
#        pass
#
#    def convert(self, position, from_position, to_position):
#        """Convert from any position type to an other.
#
#        :arg tuple position: Any position type.
#        :arg str from_position: Any position name (see _to_coordinate).
#        :arg str to_position: Any position name (see _to_position).
#
#        :return tuple: Any position type.
#        """
#        return self._to_position[to_position](
#            self._from_position[from_position](position))
