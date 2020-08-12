from .location import cut_locations, nearest_location
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

    def _nearest_region(self, coordinate):
        """
        :returns int: -1: left, 0: middle, 1: right.
        """
        outside = self._coding[1].outside(coordinate)

        if outside < 0 and self._regions[0]:
            return nearest_location(
                self._regions[:2], coordinate, self._inverted) - 1
        if outside > 0 and self._regions[2]:
            return nearest_location(
                self._regions[1:], coordinate, self._inverted)

        return 0

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

    def _coordinate_to_coding(self, coordinate):
        """Convert a coordinate to a proper coding position (c./r.).

        :arg int coordinate: Coordinate.

        :returns tuple: Coding position (c./r.).
        """
        region = self._nearest_region(coordinate)
        position = self._coding[region + 1].to_position(coordinate)

        outside = self._coding[1].orientation * region

        # Remove "outside" coordinates between regions.
        if position[2] * outside < 0:
            return (position[0], position[1], 0, outside)
        if not outside:
            _RR = self._coding[1].orientation  # FIXME: Refactor.
            if position[2] < 0 and self._regions[-1 * _RR + 1]:
                return (position[0], position[1], 0, outside)
            if position[2] > 0 and self._regions[_RR + 1]:
                return (position[0], position[1], 0, outside)

        return (*position, outside)

    def coordinate_to_coding(self, coordinate, degenerate=False):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns tuple: Coding position (c./r.).
        """
        self._check(self._coding, self._coding_error)

        position = self._coordinate_to_coding(coordinate)
        if degenerate and position[2]:
            if position[0] * position[1] < 0:
                return position[0] + position[1] - 1, 0, position[2], position[3]
            return position[0] + position[1], 0, position[2], position[3]

        return position


    def coding_to_coordinate(self, position):
        """Convert a coding position (c./r.) to a coordinate.

        :arg tuple position: Coding position (c./r.).

        :returns int: Coordinate.
        """
        self._check(self._coding, self._coding_error)

        region = self._coding[1].orientation * position[3] + 1

        # Temporary patch.
        if not self._regions[region]:
            if position[3] == 1:
                return self._coding[1].to_coordinate(
                    (position[0] + self._cds_len, position[1], 0))
            if position[3] == -1:
                return self._coding[1].to_coordinate(
                    (position[0], position[1], 0))

        return self._coding[region].to_coordinate((*position[:2], 0))

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
