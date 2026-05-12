from .multi_locus import MultiLocus


class Genomic(object):
    """Genomic crossmap object."""
    def coordinate_to_genomic(self, coordinate: int) -> dict[str, int]:
        """Convert a coordinate to a genomic position (g./m./o.).

        :arg int coordinate: Coordinate.

        :returns dict: Genomic position model.
        """
        return {'position': coordinate + 1}

    def genomic_to_coordinate(self, pos_m: dict[str, int]) -> int:
        """Convert a genomic position (g./m./o.) to a coordinate.

        :arg dict pos_m: Genomic position model.

        :returns int: Coordinate.
        """
        return pos_m['position'] - 1


class NonCoding(Genomic):
    """NonCoding crossmap object."""
    def __init__(self, locations: list[tuple[int, int]], inverted: bool = False) -> None:
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        self._inverted = inverted

        self._noncoding = MultiLocus(locations, inverted)

    def coordinate_to_noncoding(
                self,
                coordinate: int,
            ) -> dict[str, int | str]:
        """Convert a coordinate to a noncoding position (n./r.).

        :arg int coordinate: Coordinate.

        :returns dict: Noncoding position model.
        """
        multilocus_pos_m = self._noncoding.to_position(coordinate)
        return {**multilocus_pos_m, 'position': multilocus_pos_m['position'] + 1}

    def noncoding_to_coordinate(self, pos_m: dict[str, int | str]) -> int:
        """Convert a noncoding position (n./r.) to a coordinate.

        :arg dict pos_m: Noncoding position model.

        :returns int: Coordinate.
        """
        multilocus_pos_m = {**pos_m, 'position': pos_m['position'] - 1}
        return self._noncoding.to_coordinate(multilocus_pos_m)


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
        exon_end = self._noncoding.to_position(locations[-1][1] -1)

        if self._inverted:
            self._coding = (
                cds_end['position'] + cds_end['offset'],
                cds_start['position'] + cds_start['offset'] + 1
            )
            self._exons = (
                exon_end['position'] + exon_end['offset'],
                exon_start['position'] + exon_start['offset'] + 1,
            )
        else:
            self._coding = (
                cds_start['position'] + cds_start['offset'],
                cds_end['position'] + cds_end['offset'] + 1
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
        multilocus_pos_m = self._noncoding.to_position(coordinate)

        location = multilocus_pos_m['position']
        offset = multilocus_pos_m['offset']
        region = multilocus_pos_m['region']
        if region=='u':
            if self._exons[0] == self._coding[0]:
                return {
                'position': 1,
                'offset': offset,
                'region': 'u'
                }
            return {
                'position':self._coding[0],
                'offset': offset,
                'region': 'u'
            }
        if region == 'd':
            if self._exons[1] == self._coding[1]:
                return {
                'position': self._coding[1] - self._coding[0],
                'offset': offset,
                'region': 'd'
                }
            return {
                'position': location - self._coding[1] + 1,
                'offset': offset,
                'region': 'd'
            }
        if location < self._coding[0]:
            return {
                'position': self._coding[0] - location,
                'offset': offset,
                'region': '-'
            }
        if location >= self._coding[1]:
            return {
                'position': location - self._coding[1] + 1,
                'offset': offset,
                'region': '*'
            }
        return {
            'position': location - self._coding[0] + 1,
            'offset': offset,
            'region': ''
        }

    def coordinate_to_coding(self, coordinate: int, degenerate: bool = False) -> dict:
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns dict: Coding position model (c./r.).
        """
        pos_m = self._coordinate_to_coding(coordinate)

        region = pos_m['region']
        if not degenerate:
            return pos_m
        degenerate_pos_m = {**pos_m, 'offset': pos_m['offset']}
        if region == 'u':
            if self._coding[0] == 0:
                degenerate_pos_m['position'] = abs(degenerate_pos_m['offset'])
            else:
                degenerate_pos_m['position'] = degenerate_pos_m['position'] + abs(degenerate_pos_m['offset'])
            degenerate_pos_m['offset'] = 0
            degenerate_pos_m['region'] = '-'
        if region == 'd':
            if self._exons[1] == self._coding[1] :
                degenerate_pos_m['position'] = abs(degenerate_pos_m['offset'])
            else:
                degenerate_pos_m['position'] = degenerate_pos_m['position'] + abs(degenerate_pos_m['offset'])
            degenerate_pos_m['offset'] = 0
            degenerate_pos_m['region'] = '*'
        return degenerate_pos_m

    def coding_to_coordinate(self, pos_m: dict[str, int | str]) -> int:
        """Convert a coding position (c./r.) to a coordinate.

        :arg dict pos_m: Coding position model (c./r.).

        :returns int: Coordinate.
        """
        location = pos_m['position']
        region = pos_m['region']
        multilocus_pos_m = {**pos_m}

        if region in ('u', 'd'):
            multilocus_pos_m['position'] = location - 1
            return self._noncoding.to_coordinate(multilocus_pos_m)

        multilocus_pos_m['region'] = ''
        if region == '':
            multilocus_pos_m['position'] = location + self._coding[0] - 1
        elif region == '-':
            multilocus_pos_m['position'] = self._coding[0] - location
        else:
            multilocus_pos_m['position'] = self._coding[1] + location - 1
        return self._noncoding.to_coordinate(multilocus_pos_m)

    def coordinate_to_protein(self, coordinate: int) -> dict[str, int | str]:
        """Convert a coordinate to a protein position (p.).

        :arg int coordinate: Coordinate.

        :returns dict: Protein position model(p.).
        """
        pos_m = self.coordinate_to_coding(coordinate)

        location = pos_m['position']
        if pos_m['region'] in ('-', 'u'):
            return {
                'position': abs(-location // 3),
                'position_in_codon': -location % 3 + 1,
                'region': pos_m['region'],
                'offset': pos_m['offset']}
        return {
                'position': (location + 2) // 3,
                'position_in_codon': (location + 2) % 3 + 1,
                'region': pos_m['region'],
                'offset': pos_m['offset']}

    def protein_to_coordinate(self, pos_m: dict[str, int | str]) -> int:
        """Convert a protein position (p.) to a coordinate.

        :arg dict pos_m: Protein position model(p.).

        :returns int: Coordinate.
        """
        if pos_m['region'] in ('-', 'u'):
            return self.coding_to_coordinate(
                {'position': 3 * pos_m['position'] - pos_m['position_in_codon'] + 1,
                 'offset': pos_m['offset'],
                 'region': pos_m['region']})

        return self.coding_to_coordinate(
                {'position': 3 * pos_m['position'] + pos_m['position_in_codon'] - 3,
                 'offset': pos_m['offset'],
                 'region': pos_m['region']})
