from .multi_locus import MultiLocus

class Genomic(object):
    """Genomic crossmap object."""

    def coordinate_to_genomic(self, coordinate: int) -> dict[str, int]:
        """Convert a coordinate to a genomic position (g./m./o.).

        :arg int coordinate: Coordinate.

        :returns dict: Genomic position model.
        """
        return {'position': coordinate + 1}

    def genomic_to_coordinate(self, point: dict[str, int]) -> int:
        """Convert a genomic position (g./m./o.) to a coordinate.

        :arg dict point: Genomic position model.

        :returns int: Coordinate.
        """
        return point['position'] - 1


class NonCoding(Genomic):
    """NonCoding crossmap object."""

    def __init__(self, locations: list[tuple[int, int]], inverted: bool=False) -> None:
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        self._inverted = inverted

        self._noncoding = MultiLocus(locations, inverted)

    def coordinate_to_noncoding(self, coordinate: int) -> dict[str, int | str]:
        """Convert a coordinate to a noncoding position (n./r.).

        :arg int coordinate: Coordinate.

        :returns dict: Noncoding position model.
        """
        point = self._noncoding.to_position(coordinate)
        return {**point, 'position': point['position'] + 1}

    def noncoding_to_coordinate(self, point: dict[str, int | str]) -> int:
        """Convert a noncoding position (n./r.) to a coordinate.

        :arg dict point: Noncoding position model.

        :returns int: Coordinate.
        """
        return self._noncoding.to_coordinate({**point, 'position': point['position'] - 1})


class Coding(NonCoding):
    """Coding crossmap object."""
    def __init__(
            self,
            locations: list[tuple[int, int]],
            cds: tuple[int, int],
            inverted: bool = False
    ) -> None:
        """
        :arg list locations: List of locus locations.
        :arg tuple cds: Locus location.
        :arg bool inverted: Orientation.
        """
        NonCoding.__init__(self, locations, inverted)

        cds_start = self._noncoding.to_position(cds[0])
        cds_end = self._noncoding.to_position(cds[1] - 1)
        exon_start = self._noncoding.to_position(locations[0][0])
        exon_end = self._noncoding.to_position(locations[-1][1] - 1)

        if self._inverted:
            self._coding = (
                cds_end['position'] + cds_end['offset'],
                cds_start['position'] + cds_start['offset'] + 1,
            )
            self._exons = (
                exon_end['position'] + exon_end['offset'],
                exon_start['position'] + exon_start['offset'] + 1,
            )
        else:
            self._coding = (
                cds_start['position'] + cds_start['offset'],
                cds_end['position'] + cds_end['offset'] + 1,
            )
            self._exons = (
                exon_start['position'] + exon_start['offset'],
                exon_end['position'] + exon_end['offset'] + 1,
            )

    def _coordinate_to_coding(self, coordinate: int) -> dict[str, int | str]:
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.

        :returns dict: Coding position model (c./r.).
        """
        noncoding_point = self._noncoding.to_position(coordinate)

        position = noncoding_point['position']
        offset = noncoding_point['offset']
        region = noncoding_point['region']

        if region == 'u':
            if self._exons[0] == self._coding[0]:
                position = 1
            else:
                position = self._coding[0]

        elif region == 'd':
            if self._exons[1] == self._coding[1]:
                position = self._coding[1] - self._coding[0]
            else:
                position = position - self._coding[1] + 1

        elif position < self._coding[0]:
            position = self._coding[0] - position
            region = '-'

        elif position >= self._coding[1]:
            position = position - self._coding[1] + 1
            region = '*'

        else:
            position = position - self._coding[0] + 1
            region = ''

        return {
            'position': position,
            'offset': offset,
            'region': region,
        }

    def coordinate_to_coding(self, coordinate: int, degenerate: bool=False) -> dict:
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns dict: Coding position model (c./r.).
        """
        point = self._coordinate_to_coding(coordinate)

        region = point['region']
        if not degenerate:
            return point
        if region == 'u':
            if self._coding[0] == 0:
                point['position'] = abs(point['offset'])
            else:
                point['position'] = point['position'] + abs(point['offset'])
            point['offset'] = 0
            point['region'] = '-'
        elif region == 'd':
            if self._exons[1] == self._coding[1]:
                point['position'] = abs(point['offset'])
            else:
                point['position'] = point['position'] + abs(point['offset'])
            point['offset'] = 0
            point['region'] = '*'
        return point

    def coding_to_coordinate(self, point: dict[str, int | str]) -> int:
        """Convert a coding position (c./r.) to a coordinate.

        :arg dict point: Coding position model (c./r.).

        :returns int: Coordinate.
        """

        region = point['region']

        if region in ('u', 'd'):
            return self._noncoding.to_coordinate(point)

        position = point['position']
        noncoding_point = {**point, 'region': ''}
        if region == '':
            noncoding_point['position'] = position + self._coding[0] - 1
        elif region == '-':
            noncoding_point['position'] = self._coding[0] - position
        elif region == '*':
            noncoding_point['position'] = self._coding[1] + position - 1
        return self._noncoding.to_coordinate(noncoding_point)

    def coordinate_to_protein(self, coordinate: int) -> dict[str, int | str]:
        """Convert a coordinate to a protein position (p.).

        :arg int coordinate: Coordinate.

        :returns dict: Protein position model(p.).
        """
        point = self.coordinate_to_coding(coordinate)

        position = point['position']
        if point['region'] in ('-', 'u'):
            return {
                'position': abs(-position // 3),
                'position_in_codon': -position % 3 + 1,
                'region': point['region'],
                'offset': point['offset'],
            }
        return {
            'position': (position + 2) // 3,
            'position_in_codon': (position + 2) % 3 + 1,
            'region': point['region'],
            'offset': point['offset'],
        }

    def protein_to_coordinate(self, point: dict[str, int | str]) -> int:
        """Convert a protein position (p.) to a coordinate.

        :arg dict point: Protein position model(p.).

        :returns int: Coordinate.
        """
        if point['region'] in ('-', 'u'):
            return self.coding_to_coordinate(
                {
                    'position': 3 * point['position'] - point['position_in_codon'] + 1,
                    'offset': point['offset'],
                    'region': point['region'],
                }
            )

        return self.coding_to_coordinate(
            {
                'position': 3 * point['position'] + point['position_in_codon'] - 3,
                'offset': point['offset'],
                'region': point['region'],
            }
        )
