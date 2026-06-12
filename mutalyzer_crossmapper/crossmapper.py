from .multi_locus import MultiLocus
from .models import GenomicPoint, NonCodingPoint, CodingPoint,ProteinPoint


class Genomic(object):
    """Genomic crossmap object."""

    def coordinate_to_genomic(self, coordinate: int) -> dict[str, int]:
        """Convert a coordinate to a genomic point model (g./m./o.).

        :arg int coordinate: Coordinate.

        :returns dict: Genomic point model.
        """
        return GenomicPoint(coordinate + 1).to_dict()

    def genomic_to_coordinate(self, point: GenomicPoint) -> int:
        """Convert a genomic point (g./m./o.) to a coordinate.

        :arg dict point: Genomic point model.

        :returns int: Coordinate.
        """
        return GenomicPoint.to_dataclass(point).position - 1


class NonCoding(Genomic):
    """NonCoding crossmap object."""

    def __init__(self, locations: list[tuple[int, int]], inverted: bool = False) -> None:
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        self._inverted = inverted

        self._noncoding = MultiLocus(locations, inverted)

    def coordinate_to_noncoding(self, coordinate: int) -> dict[str, int | str]:
        """Convert a coordinate to a noncoding point model (n./r.).

        :arg int coordinate: Coordinate.

        :returns dict: Noncoding point model.
        """
        point = NonCodingPoint.to_dataclass(self._noncoding.to_position(coordinate))
        return NonCodingPoint(
            position=point.position + 1,
            offset=point.offset,
            region=point.region,
        ).to_dict()

    def noncoding_to_coordinate(self, point: NonCodingPoint) -> int:
        """Convert a noncoding point (n./r.) to a coordinate.

        :arg dict point: Noncoding point model.

        :returns int: Coordinate.
        """
        noncoding_point = NonCodingPoint.to_dataclass(point)
        return self._noncoding.to_coordinate(
            {
                'position': noncoding_point.position - 1,
                'region': noncoding_point.region,
                'offset': noncoding_point.offset,
            }
        )


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

    def _coordinate_to_coding(self, coordinate: int) -> CodingPoint:
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.

        :returns dict: Coding position model (c./r.).
        """
        noncoding_point = NonCodingPoint.to_dataclass(self._noncoding.to_position(coordinate))

        position = noncoding_point.position
        offset = noncoding_point.offset
        region = noncoding_point.region

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
        return CodingPoint(position=position, offset=offset, region=region)

    def coordinate_to_coding(self, coordinate: int, degenerate: bool = False) -> dict[str, int | str]:
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns dict: Coding position model (c./r.).
        """
        point = self._coordinate_to_coding(coordinate)

        region = point.region
        if not degenerate:
            return point.to_dict()

        if region == 'u':
            position = abs(point.offset) if self._coding[0] == 0 else point.position + abs(point.offset)
            return CodingPoint(position=position, offset=0, region='-').to_dict()
        if region == 'd':
            position = abs(point.offset) if self._exons[1] == self._coding[1] else point.position + abs(point.offset)
            return CodingPoint(position=position, offset=0, region='*').to_dict()
        return point.to_dict()

    def coding_to_coordinate(self, point: NonCodingPoint) -> int:
        """Convert a coding position (c./r.) to a coordinate.

        :arg dict point: Coding position model (c./r.).

        :returns int: Coordinate.
        """

        coding_point = CodingPoint.to_dataclass(point)
        region = coding_point.region

        if region in ('u', 'd'):
            return self._noncoding.to_coordinate(coding_point.to_dict())

        position = coding_point.position
        noncoding_point = {
            'position': coding_point.position,
            'region': '',
            'offset': coding_point.offset,
        }
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
        point = CodingPoint.to_dataclass(self.coordinate_to_coding(coordinate))

        position = point.position
        if point.region in ('-', 'u'):
            return ProteinPoint(
                position=abs(-position // 3),
                position_in_codon=-position % 3 + 1,
                region=point.region,
                offset=point.offset,
            ).to_dict()
        return ProteinPoint(
            position=(position + 2) // 3,
            position_in_codon=(position + 2) % 3 + 1,
            region=point.region,
            offset=point.offset,
        ).to_dict()

    def protein_to_coordinate(self, point: ProteinPoint) -> int:
        """Convert a protein position (p.) to a coordinate.

        :arg dict point: Protein position model(p.).

        :returns int: Coordinate.
        """
        protein_point = ProteinPoint.to_dataclass(point)
        if protein_point.region in ('-', 'u'):
            return self.coding_to_coordinate(
                CodingPoint(
                    position=3 * protein_point.position - protein_point.position_in_codon + 1,
                    offset=protein_point.offset,
                    region=protein_point.region,
                ).to_dict()
            )
        return self.coding_to_coordinate(
            CodingPoint(
                position=3 * protein_point.position + protein_point.position_in_codon - 3,
                offset=protein_point.offset,
                region=protein_point.region,
            ).to_dict()
        )
