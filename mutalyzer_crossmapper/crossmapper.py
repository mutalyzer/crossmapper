from .multi_locus import MultiLocus


class Genomic(object):
    """Genomic crossmap object."""
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


class NonCoding(Genomic):
    """NonCoding crossmap object."""
    def __init__(self, locations, inverted=False):
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        self._inverted = inverted

        self._noncoding = MultiLocus(locations, inverted)

    def coordinate_to_noncoding(self, coordinate):
        """Convert a coordinate to a noncoding position (n./r.).

        :arg int coordinate: Coordinate.

        :returns tuple: Noncoding position.
        """
        pos = self._noncoding.to_position(coordinate)

        return pos[0] + 1, pos[1], pos[2]

    def noncoding_to_coordinate(self, position):
        """Convert a noncoding position (n./r.) to a coordinate.

        :arg tuple position: Noncoding position.

        :returns int: Coordinate.
        """
        if position[0] > 0:
            return self._noncoding.to_coordinate(
                (position[0] - 1, position[1]))
        return self._noncoding.to_coordinate(position)


class Coding(NonCoding):
    """Coding crossmap object."""
    def __init__(self, locations, cds, inverted=False):
        """
        :arg list locations: List of locus locations.
        :arg tuple cds: Locus location.
        :arg bool inverted: Orientation.
        """
        NonCoding.__init__(self, locations, inverted)

        b0 = self._noncoding.to_position(cds[0])
        b1 = self._noncoding.to_position(cds[1])

        if self._inverted:
            self._coding = (b1[0] + b1[1] + 1, b0[0] + b0[1] + 1)
            self._cds_len = (b0[0] + b0[1]) - (b1[0] + b1[1])
        else:
            self._coding = (b0[0] + b0[1], b1[0] + b1[1])
            self._cds_len = (b1[0] + b1[1]) - (b0[0] + b0[1])

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
        pos = self._coordinate_to_coding(coordinate)

        if degenerate and pos[3]:
            if pos[2] == 0:
                if pos[0] == 1 and pos[1] < 0:
                    return pos[1], 0, -1, pos[3]
                if pos[0] == self._cds_len and pos[1] > 0:
                    return pos[0] + pos[1] - self._cds_len, 0, 1, pos[3]
            return pos[0] + pos[1], 0, pos[2], pos[3]

        return pos

    def coding_to_coordinate(self, position):
        """Convert a coding position (c./r.) to a coordinate.

        :arg tuple position: Coding position (c./r.).

        :returns int: Coordinate.
        """
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
        pos = self.coordinate_to_coding(coordinate)

        if pos[2] == -1:
            return (pos[0] // 3, pos[0] % 3 + 1, *pos[1:])
        return ((pos[0] + 2) // 3, (pos[0] + 2) % 3 + 1, *pos[1:])

    def protein_to_coordinate(self, position):
        """Convert a protein position (p.) to a coordinate.

        :arg tuple position: Protein position (p.).

        :returns int: Coordinate.
        """
        if position[3] == -1:
            return self.coding_to_coordinate(
                (3 * position[0] + position[1] - 1, *position[2:]))

        return self.coding_to_coordinate(
            (3 * position[0] + position[1] - 3, *position[2:]))
