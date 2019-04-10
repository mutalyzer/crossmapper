def _int_to_dict(i):
    """Convert an integer to a position dictionary."""
    return {'type': 'point', 'position': i}


def _tuple_to_dict(t):
    """Convert a tuple to a position dictionary."""
    if t[0] < t[4]:
        return 'upstream {} {}'.format(t[4], t[0])
    if t[0] < t[2]:
        return 'intron - {} {}'.format(t[2] + 1, t[0] - t[2])
    if t[0] < 0:
        return 'exon {}'.format(t[0])
    if t[0] > t[3]:
        return 'downstream {} {}'.format(t[3], t[0] - t[3])
    if t[0] > t[1]:
        return 'intron + {} {}'.format(t[1], t[0] - t[1])
    if t[0] > t[5]:
        return 'exon {}'.format(t[0])
    return 'cds {}'.format(t[0])


def _nearest_splice_site(coordinate, exons):
    """Given a coordinate, find the index of the nearest splice site."""
    # NOTE: For long lists, a binary search would be more efficient.
    splice_sites = [splice_site for exon in exons for splice_site in exon]

    return min(
        range(len(splice_sites)),
        key=lambda i: abs(splice_sites[i] - coordinate))


def _nearest_exon(coordinate, exons):
    return _nearest_splice_site(coordinate, exons) // 2;


def _noncoding_offsets(exons):
    upstream_length = []
    length = 0
    for exon in exons:
        upstream_length.append(length)
        length += exon[1] - exon[0]

    downstream_length = []
    length = 0
    for exon in exons[::-1]:
        downstream_length.append(length)
        length += exon[1] - exon[0]

    return list(zip(upstream_length, downstream_length[::-1]))


def _coding_offsets(exons, cds):
    upstream_length = []
    lengths = 0
    for exon in exons:
        if exon[1] >= cds[1]:
            if exon[0] <= cds[1]:
                lengths += exon[1] - cds[1]
            else:
                lengths += exon[1] - exon[0]

        upstream_length.append(lengths)

    downstream_length = []
    lengths = 0
    for exon in exons[::-1]:
        if exon[0] <= cds[0]:
            if exon[1] >= cds[0]:
                lengths += cds[0] - exon[0]
            else:
                lengths += exon[1] - exon[0]

        downstream_length.append(lengths)

    return list(zip(upstream_length, downstream_length[::-1]))


def _coordinate_to_genomic(coordinate):
    return coordinate + 1


def _genomic_to_coordinate(position):
    return position - 1


def _coordinate_to_locus(coordinate, locus, inverted=False):
    if inverted:
        if coordinate > locus[1]:
            return 1, locus[1] - coordinate
        if coordinate < locus[0]:
            return locus[1] - locus[0] + 1, locus[0] - coordinate
        return locus[1] - coordinate + 1, 0

    if coordinate < locus[0]:
        return 1, coordinate - locus[0]
    if coordinate > locus[1]:
        return locus[1] - locus[0] + 1, coordinate - locus[1]
    return coordinate - locus[0] + 1, 0


def _locus_to_coordinate(position, locus, inverted=False):
    if inverted:
        return locus[1] - position[0] - position[1] + 1
    return locus[0] + position[0] + position[1] - 1


def _coordinate_to_noncoding(coordinate, exon, offset, inverted=False):
    locus = _coordinate_to_locus(coordinate, exon, inverted)

    if inverted:
        return locus[0] + offset[1], locus[1]
    return locus[0] + offset[0], locus[1]


def _noncoding_to_coordinate(position, exon, offset, inverted=False):
    locus = _locus_to_coordinate(position, exon, inverted)

    if inverted:
        return locus - offset[1]
    return locus - offset[0]


# TODO: Position to int or tuple function.
# TODO: Transcript model to simple types function.


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


class Crossmap_(object):
    def __init__(self, transcript_model):
        """Initialise the class.

        The following assumptions must hold:

        - All positions are zero based.
        - The fields 'start', 'end' and 'type' are mandatory for all features.
        - The field 'inverted' is optional.
        - The feature of type 'cds' is optional.
        - The 'start' and 'end' positions of the transcript model are
          consistent with the first and last exon positions (if provided).
        - The 'inverted' attribute is inherited from the transcript model and
          is therefore ignored in all of its features.

        :arg dict transcript_model: Model containing exon and CDS positions.
        """
        self._inverted = transcript_model.get('inverted', False)
        self._locus = (
            transcript_model['start']['position'],
            transcript_model['end']['position'])
        self._exons = []
        self._cds = ()

        for feature in transcript_model.get('features', []):
            if feature['type'] == 'exon':
                self._exons.append(
                    (feature['start']['position'], feature['end']['position']))
            elif feature['type'] == 'cds':
                self._cds = (
                    feature['start']['position'], feature['end']['position'])

        self._to_position = {
            'genomic': self.coordinate_to_genomic,
            'coding': self.coordinate_to_coding,
            'locus': self.coordinate_to_locus,
            'noncoding': self.coordinate_to_noncoding,
            'protein': self.coordinate_to_protein}
        self._from_position = {
            'genomic': self.genomic_to_coordinate,
            'locus': self.locus_to_coordinate,
            'coding': self.coding_to_coordinate,
            'noncoding': self.noncoding_to_coordinate}

    def coordinate_to_genomic(self, coordinate):
        """Convert a zero based coordinate to a genomic position (g./m./o.).

        :arg int coordinate: Zero based coordinate.

        :returns int: Genomic position (g./m./o.).
        """
        return _coordinate_to_genomic(coordinate)

    def genomic_to_coordinate(self, position):
        """Convert a genomic position (g./m./o.) to a zero based coordinate.

        :arg int position: Genomic position (g./m./o.).

        :returns int: Coordinate.
        """
        return _genomic_to_coordinate(position)

    def coordinate_to_locus(self, coordinate):
        """Convert a zero based coordinate to a locus position.

        :arg int coordinate: Zero based coordinate.

        :returns tuple: Locus position.
        """
        return _coordinate_to_locus(coordinate, self._locus, self._inverted)

    def locus_to_coordinate(self, position):
        """Convert a locus position to a zero based coordinate.

        :arg tuple position: Locus position.

        :returns int: Coordinate.
        """
        return _locus_to_coordinate(position, self._locus, self._inverted)

    def coordinate_to_noncoding(self, coordinate):
        """Convert a zero based coordinate to a noncoding position (n./r.).

        :arg int coordinate: Zero based coordinate.

        :returns tuple: Noncoding position (n./r.).
        """
        return _coordinate_to_noncoding(
            coordinate, self._exons, self._inverted)

    def noncoding_to_coordinate(self, position):
        """Convert a noncoding position (n./r.) to a zero based coordinate.

        :arg tuple position: Noncoding position (n./r.).

        :returns int: Coordinate.
        """
        return _noncoding_to_coordinate(position, self._exons, self._inverted)

    def coordinate_to_coding(self, coordinate):
        """Convert a zero based coordinate to a coding position (c./r.).

        :arg int coordinate: Zero based coordinate.

        :returns tuple: Coding position (c./r.).
        """
        pass

    def coding_to_coordinate(self, position):
        """Convert a coding position (c./r.) to a zero based coordinate.

        :arg tuple position: Coding position (c./r.).

        :returns int: Coordinate.
        """
        pass

    def coordinate_to_protein(self, coordinate):
        """Convert a zero based coordinate to a protein position (p.).

        Note that the converse of this function does not exist.

        :arg int coordinate: Zero based coordinate.

        :returns tuple: Protein position (p.).
        """
        if not self._cds:
            raise ValueError(
                "conversion to protein position using a non coding transcript")
        pass

    def convert(self, position, from_position, to_position):
        """Convert from any position type to an other.

        :arg tuple position: Any position type.
        :arg str from_position: Any position name (see _to_coordinate).
        :arg str to_position: Any position name (see _to_position).

        :return tuple: Any position type.
        """
        return self._to_position[to_position](
            self._from_position[from_position](position))
