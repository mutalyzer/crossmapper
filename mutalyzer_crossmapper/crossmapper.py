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


def _loc(a, b):
    """Make a proper location."""
    if a >= b:
        return []
    return [(a, b)]


def _offsets(ls, inverted=False):
    """For each location, calculate the length of the preceding locations.

    :arg list ls: List of locations.
    :arg bool inverted: Direction of {ls}.

    :returns list: List of cumulative location lengths.
    """
    lengths = []

    length = 0
    direction = -1 if inverted else 1

    for location in ls[::direction]:
        lengths.append(length)
        length += location[1] - location[0]

    return lengths


def nearest_location(ls, c, p=0):
    """Find the location nearest to `c`. In case of a draw, the parameter `p`
    decides which index is chosen.

    :arg list ls: List of locations.
    :arg int c: Coordinate.
    :arg int p: Preference in case of a draw: 0: left, 1: right.
    """
    r = len(ls) - 1
    l = 0

    while l <= r:
        i = (l + r) // 2

        if c < ls[i][0]:    # `c` lies before this location.
            r = i - 1
        elif c >= ls[i][1]: # `c` lies after this location.
            l = i + 1
        else:               # `c` lies in this location.
            return i

    if c < ls[i][0]: # `c` lies before this location.
        if i and c - ls[i - 1][1] + p < ls[i][0] - c:
            return i - 1
    else:            # `c` lies after this location.
        if i < len(ls) - 1 and ls[i + 1][0] - c + p < c - ls[i][1]:
            return i + 1

    return i


def cut_locations(ls, c):
    """Divide a list of locations, cutting one of the locations in two.

    :arg int c: Coordinate.
    :arg list ls: List of locations.

    :returns tuple: locations before `c`, locations after coordinate.
    """
    i = nearest_location(ls, c)

    return ls[:i] + _loc(ls[i][0], c), _loc(c, ls[i][1]) + ls[i + 1:]


class Locus(object):
    """Locus object."""
    def __init__(self, location, inverted=False):
        """
        :arg tuple location: Locus location.
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
        index = nearest_location(self._locations, coordinate, self._inverted)
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
    """Crossmap object."""
    _noncoding_error = 'no locations provided'
    _coding_error = 'no cds provided'

    def __init__(self, locations=None, cds=None, inverted=False):
        """
        :arg list locations: List of locus locations.
        :arg tuple cds: Locus location.
        :arg bool inverted: Orientation.
        """
        self._cds = cds

        self._noncoding = None
        if locations:
            self._noncoding = MultiLocus(locations, inverted)

        self._coding = None
        self._parts = (0, 1, 2)
        if cds:
            utr5, tail = cut_locations(locations, cds[0])
            coding, utr3 = cut_locations(tail, cds[1])

            if inverted:
                utr5, utr3 = utr3, utr5
                self._parts = self._parts[::-1]

            cds_ml = MultiLocus(coding, inverted)
            self._coding = (
                MultiLocus(utr5, not inverted, True) if utr5 else cds_ml,
                cds_ml,
                MultiLocus(utr3, inverted) if utr3 else cds_ml)

    def _check(self, condition, error):
        if not condition:
            raise ValueError(error)

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

    def coordinate_to_coding(self, coordinate):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.

        :returns tuple: Coding position (c./r.).
        """
        self._check(self._coding, self._coding_error)

        if coordinate < self._cds[0]:
            return (
                *self._coding[self._parts[0]].to_position(coordinate),
                self._parts[0])
        if coordinate >= self._cds[1]:
            return (
                *self._coding[self._parts[2]].to_position(coordinate),
                self._parts[2])
        return (*self._coding[1].to_position(coordinate), 1)

    def coding_to_coordinate(self, position):
        """Convert a coding position (c./r.) to a coordinate.

        :arg tuple position: Coding position (c./r.).

        :returns int: Coordinate.
        """
        self._check(self._coding, self._coding_error)

        return self._coding[position[2]].to_coordinate(position[:2])

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
