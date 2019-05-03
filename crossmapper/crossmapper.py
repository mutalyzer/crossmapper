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


def _nearest_exon(coordinate, boundaries):
    return _nearest_boundary(coordinate, boundaries) // 2;


# TODO: From position and offset to nearest exon.


def _cut(coordinate, exons):
    """Divide a list of exons, cutting one of the exons in two.

    :arg int coordinate: Coordinate.
    :arg list exons: List of exons.

    :returns tuple: {exons} before coordinate, {exons} after coordinate.
    """
    exon = _nearest_exon(coordinate, exons)

    return (
        exons[:exon] + [(exons[exon][0], coordinate)],
        [(coordinate, exons[exon][1])] + exons[exon + 1:])


def _negate(position):
    return (-position[0], -position[1])


def _offsets(exons, inverted=False):
    """For each exon, calculate the length of the preceding exons.

    :arg list exons: List of exons.
    :arg bool inverted: Direction of {exons}.

    :returns list: List of cumulative exon lengths.
    """
    direction = 1
    if inverted:
        direction = -1

    length = 0
    lengths = []
    for exon in exons[::direction]:
        lengths.append(length)
        length += exon[1] - exon[0]

    return lengths[::direction]


def _coordinate_to_genomic(coordinate):
    return coordinate + 1


def _genomic_to_coordinate(position):
    return position - 1


def _coordinate_to_locus(coordinate, locus, inverted=False):
    """Convert a coordinate to a position relative to a locus.

    :arg int coordinate: Coordinate.
    :arg tuple locus: Locus.
    :arg bool inverted: Direction of {locus}.

    :returns tuple: Position relative to {locus}.
    """
    if inverted:
        if coordinate >= locus[1]:
            return 1, locus[1] - coordinate - 1
        if coordinate < locus[0]:
            return locus[1] - locus[0], locus[0] - coordinate
        return locus[1] - coordinate, 0

    if coordinate < locus[0]:
        return 1, coordinate - locus[0]
    if coordinate >= locus[1]:
        return locus[1] - locus[0], coordinate - locus[1] + 1
    return coordinate - locus[0] + 1, 0


def _locus_to_coordinate(position, locus, inverted=False):
    """Convert a position relative to a locus to a coordinate.

    :arg int position: Position.
    :arg tuple locus: Locus.
    :arg bool inverted: Direction of {locus}.

    :returns int: Coordinate.
    """
    if inverted:
        return locus[1] - position[0] - position[1]
    return locus[0] + position[0] + position[1] - 1


def _coordinate_to_exon(coordinate, exon, offset, inverted=False):
    """Convert a coordinate to a position relative to an exon.

    :arg int coordinate: Coordinate.
    :arg tuple exon: Exon.
    :arg int offset: Offset of {exon}.
    :arg bool inverted: Direction of {exon}.

    :returns tuple: Position relative to {exon}.
    """
    locus = _coordinate_to_locus(coordinate, exon, inverted)

    return locus[0] + offset, locus[1]


def _exon_to_coordinate(position, exon, offset, inverted=False):
    """Convert a position relative to an exon to a coordinate.

    :arg int position: Position.
    :arg tuple exon: Exon.
    :arg int offset: Offset of {exon}.
    :arg bool inverted: Direction of {exon}.

    :returns int: Coordinate.
    """
    return _locus_to_coordinate(
        (position[0] - offset, position[1]), exon, inverted)


class Crossmap(object):
    def __init__(self, locus, exons=None, cds=None, inverted=False):
        self._inverted = inverted
        self._locus = locus
        self._exons = exons
        self._cds = cds
        self._inverted = inverted

        self._boundaries = None
        self._offsets = None
        if self._exons:
            self._boundaries = sum([e[0], e[1] - 1] for e in self._exons, [])
            self._offsets = _offsets(self._exons)

        self._cds_exons = None
        self._cds_offsets = None
        if self._cds:
            utr5, tail = cut(self._cds[0], self._exons)
            coding, utr3 = cut(self._cds[1], tail)
            self._cds_exons = (utr5, coding, utr3)
            self._cds_offsets = (
                _offsets(utr5, True), _offsets(coding), _offsets(utr3))

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
        exon = _nearest_exon(coordinate, self._boundaries)

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
