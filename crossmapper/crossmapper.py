class Crossmap(object):
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

        :arg dict coordinate: Zero based coordinate.

        :returns dict: Genomic position (g./m./o.).
        """
        return {'position': coordinate['position'] + 1}

    def genomic_to_coordinate(self, genomic_position):
        """Convert a genomic position (g./m./o.) to a zero based coordinate.

        :arg dict genomic_position: Genomic position (g./m./o.).

        :returns dict: Coordinate.
        """
        return {'position': genomic_position['position'] - 1}

    def coordinate_to_locus(self, coordinate):
        """Convert a zero based coordinate to a locus position.

        :arg dict coordinate: Zero based coordinate.

        :returns dict: Locus position.
        """
        if self._inverted:
            return {
                'position': self._locus[1] + 1,
                'offset': {
                    'position': self._locus[1] - coordinate['position']}}
        return {
            'position': self._locus[0] + 1,
            'offset': {
                'position': coordinate['position'] - self._locus[0]}}

    def locus_to_coordinate(self, locus_position):
        """Convert a locus position to a zero based coordinate.

        :arg dict locus_position: Locus position.

        :returns dict: Coordinate.
        """
        if self._inverted:
            return {
                'position': locus_position['position'] -
                    locus_position['offset']['position'] - 1}
        return {
            'position': locus_position['position'] +
                locus_position['offset']['position'] - 1}

    def coordinate_to_noncoding(self, coordinate):
        """Convert a zero based coordinate to a noncoding position (n./r.).

        :arg dict coordinate: Zero based coordinate.

        :returns dict: Noncoding position (n./r.).
        """
        pass

    def noncoding_to_coordinate(self, noncoding_position):
        """Convert a noncoding position (n./r.) to a zero based coordinate.

        :arg dict noncoding_position: Noncoding position (n./r.).

        :returns dict: Coordinate.
        """
        pass

    def coordinate_to_coding(self, coordinate):
        """Convert a zero based coordinate to a coding position (c./r.).

        :arg dict coordinate: Zero based coordinate.

        :returns dict: Coding position (c./r.).
        """
        pass

    def coding_to_coordinate(self, coding_position):
        """Convert a coding position (c./r.) to a zero based coordinate.

        :arg dict coding_position: Coding position (c./r.).

        :returns dict: Coordinate.
        """
        pass

    def coordinate_to_protein(self, coordinate):
        """Convert a zero based coordinate to a protein position (p.).

        Note that the converse of this function does not exist.

        :arg dict coordinate: Zero based coordinate.

        :returns dict: Protein position (p.).
        """
        if not self._cds:
            raise ValueError(
                "conversion to protein position using a non coding transcript")
        pass

    def convert(self, position, from_position, to_position):
        """Convert from any position type to an other.

        :arg dict position: Any position type.
        :arg str from_position: Any position name (see _to_coordinate).
        :arg str to_position: Any position name (see _to_position).

        :return dict: Any position type.
        """
        return self._to_position[to_position](
            self._from_position[from_position](position))
