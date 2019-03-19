class Crossmap(object):
    def __init__(self, transcript_model):
        """Initialise the class.

        The following assumptions must hold:

        - All positions are zero based.
        - The fields 'start', 'end' and 'type' are mandatory for all features.
        - The field 'inverted' is optional.
        - The 'start' and 'end' positions of the transcript model are
          consistent with the first and last exon positions and are therefore
          ignored.
        - The 'inverted' attribute is inherited from the transcript model and
          is therefore ignored in all of its features.

        :arg dict transcript_model: Model containing exon and cds positions.
        """
        self._inverted = transcript_model.get('inverted', False)
        self._exons = []
        self._cds = ()

        for feature in transcript_model['features']:
            if feature['type'] == 'exon':
                self._exons.append(
                    (feature['start']['position'], feature['end']['position']))
            elif feature['type'] == 'cds':
                self._cds = (
                    feature['start']['position'], feature['end']['position'])

        self._to_position = {
            'genomic': self.coordinate_to_genomic,
            'coding': self.coordinate_to_coding,
            'noncoding': self.coordinate_to_noncoding,
            'protein': self.coordinate_to_protein}
        self._from_position = {
            'genomic': self.genomic_to_coordinate,
            'coding': self.coding_to_coordinate,
            'noncoding': self.noncoding_to_coordinate}

    def coordinate_to_genomic(self, coordinate):
        """Convert a zero based coordinate to a genomic position (g./m./o.).

        :arg dict coordinate: Zero based coordinate.

        :returns dict: Genomic position (g./m./o.).
        """
        return {'type': 'point', 'position': coordinate['position'] + 1}

    def genomic_to_coordinate(self, genomic_position):
        """Convert a genomic position (g./m./o.) to a zero based coordinate.

        :arg dict genomic_position: Genomic position (g./m./o.).

        :returns dict: Coordinate.
        """
        return {'type': 'point', 'position': genomic_position['position'] - 1}

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

    def coordinate_to_protein(self, coordinate):
        """Convert a zero based coordinate to a protein position (p.).

        Note that the converse of this function does not exist.

        :arg dict coordinate: Zero based coordinate.

        :returns dict: Protein position (p.).
        """
        pass

    def convert(self, position, from_position, to_position):
        """Convert from any position type to an other.

        :arg dict position: Any position.
        :arg str from_position: Any position name (see _to_coordinate).
        :arg str to_position: Any position name (see _to_position).

        :return dict: Any position.
        """
        return self._to_position[to_position](
            self._from_position[from_position](position))


transcript_model = {
    'features': [
        {
            'start': {
                'position': 10},
            'end': {
                'position': 20},
            'type': 'exon'},
        {
            'start': {
                'position': 30},
            'end': {
                'position': 40},
            'type': 'exon'},
        {
            'start': {
                'position': 50},
            'end': {
                'position': 60},
            'type': 'exon'},
        {
            'start': {
                'position': 70},
            'end': {
                'position': 80},
            'type': 'exon'},
        {
            'start': {
                'position': 35},
            'end': {
                'position': 75},
            'type': 'cds'}],
    'start': {
        'position': 10},
    'end': {
        'position': 75},
    'inverted': True,
    'type': 'rna'}

location = {
    'type': 'point',
    'position': 100,
    'offset': {
        'value': 3}}

crossmap = Crossmap(transcript_model)
print(crossmap._inverted)
print(crossmap._exons)
print(crossmap._cds)

print(crossmap.coordinate_to_genomic(location))
print(crossmap.genomic_to_coordinate(location))

print(crossmap.convert(location, 'genomic', 'genomic'))
