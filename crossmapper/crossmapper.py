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


def _nearest_locus(coordinate, boundaries):
    return _nearest_boundary(coordinate, boundaries) // 2


def _cut(coordinate, loci):
    """Divide a list of loci, cutting one of the loci in two.

    :arg int coordinate: Coordinate.
    :arg list loci: List of loci.

    :returns tuple: {loci} before coordinate, {loci} after coordinate.
    """
    locus = _nearest_locus(coordinate, loci)

    return (
        loci[:locus] + [(loci[locus][0], coordinate)],
        [(coordinate, loci[locus][1])] + loci[locus + 1:])


#def _negative(position):
#    return (-position[0], -position[1])


def _offsets(loci, inverted=False):
    """For each locus, calculate the length of the preceding loci.

    :arg list loci: List of loci.
    :arg bool inverted: Direction of {loci}.

    :returns list: List of cumulative locus lengths.
    """
    direction = 1
    if inverted:
        direction = -1

    length = 0
    lengths = []
    for locus in loci[::direction]:
        lengths.append(length)
        length += locus[1] - locus[0]

    return lengths[::direction]


def to_position(coordinate):
    return coordinate + 1


def to_coordinate(position):
    return position - 1


#def _coordinate_to_locus(coordinate, locus, inverted=False):
#    """Convert a coordinate to a position relative to a locus.
#
#    :arg int coordinate: Coordinate.
#    :arg tuple locus: Locus.
#    :arg bool inverted: Direction of {locus}.
#
#    :returns tuple: Position relative to {locus}.
#    """
#    if inverted:
#        if coordinate >= locus[1]:
#            return 1, locus[1] - coordinate - 1
#        if coordinate < locus[0]:
#            return locus[1] - locus[0], locus[0] - coordinate
#        return locus[1] - coordinate, 0
#
#    if coordinate < locus[0]:
#        return 1, coordinate - locus[0]
#    if coordinate >= locus[1]:
#        return locus[1] - locus[0], coordinate - locus[1] + 1
#    return coordinate - locus[0] + 1, 0


#def _locus_to_coordinate(position, locus, inverted=False):
#    """Convert a position relative to a locus to a coordinate.
#
#    :arg int position: Position.
#    :arg tuple locus: Locus.
#    :arg bool inverted: Direction of {locus}.
#
#    :returns int: Coordinate.
#    """
#    if inverted:
#        return locus[1] - position[0] - position[1]
#    return locus[0] + position[0] + position[1] - 1


#def _coordinate_to_exon(coordinate, exon, offset, inverted=False):
#    """Convert a coordinate to a position relative to an exon.
#
#    :arg int coordinate: Coordinate.
#    :arg tuple exon: Exon.
#    :arg int offset: Offset of {exon}.
#    :arg bool inverted: Direction of {exon}.
#
#    :returns tuple: Position relative to {exon}.
#    """
#    locus = _coordinate_to_locus(coordinate, exon, inverted)
#
#    return locus[0] + offset, locus[1]


#def _exon_to_coordinate(position, exon, offset, inverted=False):
#    """Convert a position relative to an exon to a coordinate.
#
#    :arg int position: Position.
#    :arg tuple exon: Exon.
#    :arg int offset: Offset of {exon}.
#    :arg bool inverted: Direction of {exon}.
#
#    :returns int: Coordinate.
#    """
#    return _locus_to_coordinate(
#        (position[0] - offset, position[1]), exon, inverted)


class Locus(object):
    def __init__(self, locus, inverted=False):
        """Construct a Locus object.

        :arg tuple locus: Locus coordinates.
        :arg bool inverted: Orientation of {locus}.
        """
        self._locus = locus
        self._inverted = inverted

    def to_position(self, coordinate):
        """Convert a coordinate to a position.

        :arg int coordinate: Coordinate.

        :returns tuple: Position.
        """
        if self._inverted:
            if coordinate >= self._locus[1]:
                return 1, self._locus[1] - coordinate - 1
            if coordinate < self._locus[0]:
                return (
                    self._locus[1] - self._locus[0],
                    self._locus[0] - coordinate)
            return self._locus[1] - coordinate, 0

        if coordinate < self._locus[0]:
            return 1, coordinate - self._locus[0]
        if coordinate >= self._locus[1]:
            return (
                self._locus[1] - self._locus[0],
                coordinate - self._locus[1] + 1)
        return coordinate - self._locus[0] + 1, 0

    def to_coordinate(self, position):
        """Convert a position to a coordinate.

        :arg int position: Position.

        :returns int: Coordinate.
        """
        if self._inverted:
            return self._locus[1] - position[0] - position[1]
        return self._locus[0] + position[0] + position[1] - 1


class MultiLocus(object):
    def __init__(self, locus_list, inverted=False, negated=False):
        """Construct a MultiLocus object.
        :arg list locus_list: List of locus coordinates.
        :arg bool inverted: Orientation of {locus}.
        :arg bool negated: Change the sign of all positions.
        """
        self._inverted = inverted
        self._negated = negated

        self._loci = [Locus(locus, inverted) for locus in locus_list]
        self._boundaries = sum(
            [[locus[0], locus[1] - 1] for locus in locus_list], [])

        self._offsets = _offsets(locus_list, inverted)
        if inverted:
            self._offsets = self._offsets[::-1]

    def _sign(self, position):
        if self._negated:
            return -position[0], -position[1]
        return position

    def to_position(self, coordinate):
        """Convert a coordinate to a position.

        :arg int coordinate: Coordinate.

        :returns tuple: Position.
        """
        index = _nearest_boundary(coordinate, self._boundaries) // 2
        locus = self._loci[index].to_position(coordinate)

        return self._sign((locus[0] + self._offsets[index], locus[1]))

    def to_coordinate(self, position):
        """Convert a position to a coordinate.

        :arg int position: Position.

        :returns int: Coordinate.
        """
        position_ = self._sign(position)

        index = min(
            len(self._offsets),
            max(0, bisect_left(self._offsets, position_[0]) - 1))

        return self._loci[index].to_coordinate(
            (position_[0] - self._offsets[index], position_[1]))


class Crossmap(object):
    def __init__(self, locus=None, exons=None, cds=None, inverted=False):
        self._locus = None
        self._noncoding = None
        self._coding = None

        if locus:
            self._locus = Locus(locus, inverted)
        if exons:
            self._noncoding = MultiLocus(exons, inverted)
        if cds:
            utr5, tail = cut(cds[0], self._exons)
            coding, utr3 = cut(cds[1], tail)

            self._coding = (
                MultiLocus(utr5, inverted, True),
                MultiLocus(coding, inverted),
                MultiLocus(utr3, inverted))

    def coordinate_to_genomic(self, coordinate):
        """Convert a coordinate to a genomic position (g./m./o.).

        :arg int coordinate: Coordinate.

        :returns int: Genomic position.
        """
        return _coordinate_to_genomic(coordinate)

    def genomic_to_coordinate(self, position):
        """Convert a genomic position (g./m./o.) to a coordinate.

        :arg int position: Genomic position.

        :returns int: Coordinate.
        """
        return _genomic_to_coordinate(position)

    def coordinate_to_locus(self, coordinate):
        """Convert a coordinate to a locus position.

        :arg int coordinate: Coordinate.

        :returns tuple: Locus position.
        """
        return _coordinate_to_locus(coordinate, self._locus, self._inverted)

    def locus_to_coordinate(self, position):
        """Convert a locus position to a coordinate.

        :arg tuple position: Locus position.

        :returns int: Coordinate.
        """
        return _locus_to_coordinate(position, self._locus, self._inverted)

    def coordinate_to_noncoding(self, coordinate):
        """Convert a coordinate to a noncoding position (n./r.).

        :arg int coordinate: Coordinate.

        :returns tuple: Noncoding position.
        """
        exon = _nearest_locus(coordinate, self._boundaries)

        return _coordinate_to_exon(
            coordinate, self._exons[exon], self._offsets[exon], self._inverted)

    #def noncoding_to_coordinate(self, position):
    #    """Convert a noncoding position (n./r.) to a coordinate.

    #    :arg tuple position: Noncoding position.

    #    :returns int: Coordinate.
    #    """
    #    return _exon_to_coordinate(position, self._exons, self._inverted)

#    def coordinate_to_coding(self, coordinate):
#        """Convert a coordinate to a coding position (c./r.).
#
#        :arg int coordinate: Coordinate.
#
#        :returns tuple: Coding position (c./r.).
#        """
#        pass
#
#    def coding_to_coordinate(self, position):
#        """Convert a coding position (c./r.) to a coordinate.
#
#        :arg tuple position: Coding position (c./r.).
#
#        :returns int: Coordinate.
#        """
#        pass
#
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
