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


def _nearest_boundary(coordinate, exons):
    """Given a coordinate, find the index of the nearest exon boundary.

    When a draw occurs, the right boundary is chosen over the left one.

    :arg int coordinate: Coordinate.
    :arg list exons: List of exons.

    :returns int: Index of nearest boundary.
    """
    # NOTE: A binary search would be more efficient.
    # FIXME: Breaks down for adjacent exons.
    boundaries = [boundary for exon in exons for boundary in exon]

    return min(
        range(len(boundaries)),
        key=lambda i: abs(boundaries[i] - coordinate - 1))


def _nearest_exon(coordinate, exons):
    return _nearest_boundary(coordinate, exons) // 2;


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

    :arg int coordinate: Coordinate.
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

    :arg int coordinate: Coordinate.
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

        self._cummulative_exon_lengtht = None # Name them offsets.
        if (self._exons):
            self._cummulative_exon_lengtht = _cummulative_exon_lengtht(
                self._exons)

        self._cummulative_non_cds_lengtht = None
        if (self._cds):
            self._cummulative_non_cds_lengtht = _cummulative_non_cds_lengtht(
                self._exons, self._cds)


#class Crossmap_(object):
#    def __init__(self, transcript_model):
#        """Initialise the class.
#
#        The following assumptions must hold:
#
#        - All positions are zero based.
#        - The fields 'start', 'end' and 'type' are mandatory for all features.
#        - The field 'inverted' is optional.
#        - The feature of type 'cds' is optional.
#        - The 'start' and 'end' positions of the transcript model are
#          consistent with the first and last exon positions (if provided).
#        - The 'inverted' attribute is inherited from the transcript model and
#          is therefore ignored in all of its features.
#
#        :arg dict transcript_model: Model containing exon and CDS positions.
#        """
#        self._inverted = transcript_model.get('inverted', False)
#        self._locus = (
#            transcript_model['start']['position'],
#            transcript_model['end']['position'])
#        self._exons = []
#        self._cds = ()
#
#        for feature in transcript_model.get('features', []):
#            if feature['type'] == 'exon':
#                self._exons.append(
#                    (feature['start']['position'], feature['end']['position']))
#            elif feature['type'] == 'cds':
#                self._cds = (
#                    feature['start']['position'], feature['end']['position'])
#
#        self._to_position = {
#            'genomic': self.coordinate_to_genomic,
#            'coding': self.coordinate_to_coding,
#            'locus': self.coordinate_to_locus,
#            'noncoding': self.coordinate_to_noncoding,
#            'protein': self.coordinate_to_protein}
#        self._from_position = {
#            'genomic': self.genomic_to_coordinate,
#            'locus': self.locus_to_coordinate,
#            'coding': self.coding_to_coordinate,
#            'noncoding': self.noncoding_to_coordinate}
#
#    def coordinate_to_genomic(self, coordinate):
#        """Convert a zero based coordinate to a genomic position (g./m./o.).
#
#        :arg int coordinate: Zero based coordinate.
#
#        :returns int: Genomic position (g./m./o.).
#        """
#        return _coordinate_to_genomic(coordinate)
#
#    def genomic_to_coordinate(self, position):
#        """Convert a genomic position (g./m./o.) to a zero based coordinate.
#
#        :arg int position: Genomic position (g./m./o.).
#
#        :returns int: Coordinate.
#        """
#        return _genomic_to_coordinate(position)
#
#    def coordinate_to_locus(self, coordinate):
#        """Convert a zero based coordinate to a locus position.
#
#        :arg int coordinate: Zero based coordinate.
#
#        :returns tuple: Locus position.
#        """
#        return _coordinate_to_locus(coordinate, self._locus, self._inverted)
#
#    def locus_to_coordinate(self, position):
#        """Convert a locus position to a zero based coordinate.
#
#        :arg tuple position: Locus position.
#
#        :returns int: Coordinate.
#        """
#        return _locus_to_coordinate(position, self._locus, self._inverted)
#
#    def coordinate_to_noncoding(self, coordinate):
#        """Convert a zero based coordinate to a noncoding position (n./r.).
#
#        :arg int coordinate: Zero based coordinate.
#
#        :returns tuple: Noncoding position (n./r.).
#        """
#        return _coordinate_to_noncoding(
#            coordinate, self._exons, self._inverted)
#
#    def noncoding_to_coordinate(self, position):
#        """Convert a noncoding position (n./r.) to a zero based coordinate.
#
#        :arg tuple position: Noncoding position (n./r.).
#
#        :returns int: Coordinate.
#        """
#        return _noncoding_to_coordinate(position, self._exons, self._inverted)
#
#    def coordinate_to_coding(self, coordinate):
#        """Convert a zero based coordinate to a coding position (c./r.).
#
#        :arg int coordinate: Zero based coordinate.
#
#        :returns tuple: Coding position (c./r.).
#        """
#        pass
#
#    def coding_to_coordinate(self, position):
#        """Convert a coding position (c./r.) to a zero based coordinate.
#
#        :arg tuple position: Coding position (c./r.).
#
#        :returns int: Coordinate.
#        """
#        pass
#
#    def coordinate_to_protein(self, coordinate):
#        """Convert a zero based coordinate to a protein position (p.).
#
#        Note that the converse of this function does not exist.
#
#        :arg int coordinate: Zero based coordinate.
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
