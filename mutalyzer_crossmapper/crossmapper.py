from .multi_locus import MultiLocus


class Crossmap(object):
    """Crossmap object."""
    _noncoding_error = 'no locations provided'
    _coding_error = 'no cds provided'
    _position_error = 'invalid position'

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
                b0 = self._noncoding.to_position(cds[0])
                b1 = self._noncoding.to_position(cds[1])

                if self._inverted:
                    self._coding = (b1[0] + b1[1] + 1, b0[0] + b0[1] + 1)
                else:
                    self._coding = (b0[0] + b0[1], b1[0] + b1[1])

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

        pos = self._noncoding.to_position(coordinate)

        if pos[0] >= 0:
            return pos[0] + 1, pos[1], pos[2]
        return pos

    def noncoding_to_coordinate(self, position):
        """Convert a noncoding position (n./r.) to a coordinate.

        :arg tuple position: Noncoding position.

        :returns int: Coordinate.
        """
        self._check(self._noncoding, self._noncoding_error)
        self._check(position[0], self._position_error)

        if position[0] > 0:
            return self._noncoding.to_coordinate(
                (position[0] - 1, position[1]))
        return self._noncoding.to_coordinate(position)

    def _coordinate_to_coding(self, coordinate):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.

        :returns tuple: Coding position (c./r.).
        """
        pos = self._noncoding.to_position(coordinate)

        if pos[0] < self._coding[0]:
            return pos[0] - self._coding[0], pos[1], -1, pos[2]
        elif pos[0] >= self._coding[1]:
            return pos[0] - self._coding[1] + 1, pos[1], 1, pos[2]
        return pos[0] - self._coding[0] + 1, pos[1], 0, pos[2]

    def coordinate_to_coding(self, coordinate, degenerate=False):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns tuple: Coding position (c./r.).
        """
        self._check(self._coding, self._coding_error)

        pos = self._coordinate_to_coding(coordinate)

        if degenerate and pos[3]:
            region = 1 if pos[3] > 0 else -1

            if pos[0] == 1 and pos[1] < 0:
                return pos[1], 0, region, pos[3]
            return pos[0] + pos[1], 0, region, pos[3]

        return pos

    def coding_to_coordinate(self, position):
        """Convert a coding position (c./r.) to a coordinate.

        :arg tuple position: Coding position (c./r.).

        :returns int: Coordinate.
        """
        self._check(self._coding, self._coding_error)
        self._check(position[0], self._position_error)

        if position[2] == -1:
            return self._noncoding.to_coordinate(
                (position[0] + self._coding[0], position[1]))
        elif position[2] == 1:
            return self._noncoding.to_coordinate(
                (position[0] + self._coding[1] - 1, position[1]))
        return self._noncoding.to_coordinate(
            (position[0] + self._coding[0] - 1, position[1]))

    def coordinate_to_protein(self, coordinate):
        """Convert a coordinate to a protein position (p.).

        :arg int coordinate: Coordinate.

        :returns tuple: Protein position (p.).
        """
        self._check(self._coding, self._coding_error)

        pos = self.coordinate_to_coding(coordinate)

        if not pos[2]:
            return pos[0] // 3, pos[0] % 3 + 1, pos[1], 0
        return (pos[0] + 2) // 3, (pos[0] + 2) % 3 + 1, pos[1], pos[2]

    def protein_to_coordinate(self, position):
        """Convert a protein position (p.) to a coordinate.

        :arg tuple position: Protein position (p.).

        :returns int: Coordinate.
        """
        self._check(self._coding, self._coding_error)
        self._check(position[0], self._position_error)

        if not position[3]:
            return self.coding_to_coordinate(
                (position[0] * 3 + position[1] - 1, position[2], position[3]))

        return self.coding_to_coordinate(
            (position[0] * 3 - 3 + position[1], position[2], position[3]))
