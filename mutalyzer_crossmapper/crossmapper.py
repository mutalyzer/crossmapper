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


def _nearest_boundary(lb, rb, c, p):
    """Find the boundary nearest to `c`. In case of a draw, the parameter `p`
    decides which one is chosen.

    :arg int lb: Left boundary.
    :arg int rb: Right boundary.
    :arg int c: Coordinate (`lb` <= `c` <= `rb`)).
    :arg int p: Preference in case of a draw: 0: left, 1: right.

    :returns int: Nearest boundary: 0: left, 1: right.
    """
    dl = c - lb + 1
    dr = rb - c

    if dl < dr:
        return 0
    if dl > dr:
        return 1
    return p


def nearest_location(ls, c, p=0):
    """Find the location nearest to `c`. In case of a draw, the parameter `p`
    decides which index is chosen.

    :arg list ls: List of locations.
    :arg int c: Coordinate.
    :arg int p: Preference in case of a draw: 0: left, 1: right.

    :returns int: Nearest location.
    """
    rb = len(ls) - 1
    lb = 0

    while lb <= rb:
        i = (lb + rb) // 2

        if c < ls[i][0]:     # `c` lies before this location.
            rb = i - 1
        elif c >= ls[i][1]:  # `c` lies after this location.
            lb = i + 1
        else:                # `c` lies in this location.
            return i

    if i and c < ls[i][0]:  # `c` lies before this location.
        return i - 1 + _nearest_boundary(ls[i - 1][1], ls[i][0], c, p)
    if i < len(ls) - 1:     # `c` lies after this location.
        return i + _nearest_boundary(ls[i][1], ls[i + 1][0], c, p)

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
            if position[0] < 1:  # Degenerate position.
                return self._location[1] - position[0] - position[1] - 1
            return self._location[1] - position[0] - position[1]

        if position[0] < 1:      # Degenerate position.
            return self._location[0] + position[0] + position[1]
        return self._location[0] + position[0] + position[1] - 1


class MultiLocus(object):
    """MultiLocus object."""
    # TODO: Add region coordinate, also to Crossmap.
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

    def _direction(self, index):
        if self._inverted:
            return 2 - index
        return index

    def _nearest_region(self, coordinate):
        if coordinate < self._regions[1][0]:
            if self._regions[0]:
                return nearest_location(
                    self._regions[:2], coordinate, self._inverted)
        if coordinate >= self._regions[1][1]:
            if self._regions[2]:
                return nearest_location(
                    self._regions[1:], coordinate, self._inverted) + 1
        return 1

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
        position = self._coding[region].to_position(coordinate)
        selected_region = self._direction(region)

        if (
                degenerate and
                not self._regions[1][0] <= coordinate < self._regions[1][1]):
            if selected_region == 1:
                if (
                        not self._regions[self._direction(0)] and
                        position[0] == 1 and position[1] < 0):
                    return (position[1], 0, 0)
                if (
                        not self._regions[self._direction(2)] and
                        position[0] == self._cds_len and position[1] > 0):
                    return (position[1], 0, 2)
            return (position[0] + position[1], 0, selected_region)

        return (*position, selected_region)

    def coding_to_coordinate(self, position):
        """Convert a coding position (c./r.) to a coordinate.

        :arg tuple position: Coding position (c./r.).

        :returns int: Coordinate.
        """
        self._check(self._coding, self._coding_error)

        region = self._direction(position[2])
        if not self._regions[region]:  # Degenerate position.
            if position[2] == 0:
                return self._coding[1].to_coordinate(position[:2])
            if position[2] == 2:
                return self._coding[1].to_coordinate(
                    (self._cds_len + position[0], position[1], 1))

        return self._coding[region].to_coordinate(position[:2])

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
