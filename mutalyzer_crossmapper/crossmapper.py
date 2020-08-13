from .multi_locus import MultiLocus


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

        self._noncoding = ()
        self._coding = ()
        if locations:
            self._noncoding = MultiLocus(locations, inverted)

            if cds:
                self._coding = (
                    self._noncoding.to_position(cds[0])[0] + self._noncoding.to_position(cds[0])[1],
                    self._noncoding.to_position(cds[1])[0] + self._noncoding.to_position(cds[1])[1])
                if self._inverted:
                    self._coding = self._coding[1] + 1, self._coding[0] + 1

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

    def coordinate_to_coding(self, coordinate, degenerate=False):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns tuple: Coding position (c./r.).
        """
        self._check(self._coding, self._coding_error)

        position = self._noncoding.to_position(coordinate, degenerate)

        if position[0] < self._coding[0]:
            return position[0] - self._coding[0], position[1], position[2], -1
        elif position[0] >= self._coding[1]:
            return position[0] - self._coding[1] + 1, position[1], position[2], 1
        return position[0] - self._coding[0] + 1, position[1], position[2], 0

    def coding_to_coordinate(self, position):
        """Convert a coding position (c./r.) to a coordinate.

        :arg tuple position: Coding position (c./r.).

        :returns int: Coordinate.
        """
        self._check(self._coding, self._coding_error)

        if position[3] == -1:
            return self._noncoding.to_coordinate(
                (position[0] + self._coding[0], *position[1:]))
        elif position[3] == 1:
            return self._noncoding.to_coordinate(
                (position[0] + self._coding[1] - 1, *position[1:]))
        return self._noncoding.to_coordinate(
            (position[0] + self._coding[0] - 1, *position[1:]))

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
