from .multi_locus import MultiLocus
from .models import GenomicPoint, NonCodingPoint, CodingPoint, ProteinPoint


class Genomic(object):
    """Genomic crossmap object."""

    def coordinate_to_genomic(self, coordinate: int) -> GenomicPoint:
        """Convert a coordinate to a genomic point model (g./m./o.).

        :arg int coordinate: Coordinate.

        :returns GenomicPoint: Genomic point model.
        """
        return GenomicPoint(coordinate + 1)

    def genomic_to_coordinate(self, point: GenomicPoint) -> int:
        """Convert a genomic point (g./m./o.) to a coordinate.

        :arg GenomicPoint point: Genomic point model.

        :returns int: Coordinate.
        """
        return point.position - 1


class NonCoding(Genomic):
    """NonCoding crossmap object."""

    def __init__(self, locations: list[tuple[int, int]], inverted: bool = False) -> None:
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        self._inverted = inverted

        self._noncoding = MultiLocus(locations, inverted)

    def coordinate_to_noncoding(self, coordinate: int) -> NonCodingPoint:
        """Convert a coordinate to a noncoding point model (n./r.).

        :arg int coordinate: Coordinate.

        :returns NonCodingPoint: Noncoding point model.
        """
        point = self._noncoding.to_position(coordinate)
        return NonCodingPoint(
            position=point.position + 1,
            offset=point.offset,
            region=point.region
        )

    def noncoding_to_coordinate(self, point: NonCodingPoint) -> int:
        """Convert a noncoding point (n./r.) to a coordinate.

        :arg NonCodingPoint point: Noncoding point model.

        :returns int: Coordinate.
        """
        return self._noncoding.to_coordinate(
            NonCodingPoint(
                position=point.position - 1,
                offset=point.offset,
                region=point.region
            )
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
                cds_end.position + cds_end.offset,
                cds_start.position + cds_start.offset + 1
            )
            self._exons = (
                exon_end.position + exon_end.offset,
                exon_start.position + exon_start.offset + 1
            )
        else:
            self._coding = (
                cds_start.position + cds_start.offset,
                cds_end.position + cds_end.offset + 1
            )
            self._exons = (
                exon_start.position + exon_start.offset,
                exon_end.position + exon_end.offset + 1
            )

    def _coordinate_to_coding(self, coordinate: int) -> CodingPoint:
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.

        :returns CodingPoint: Coding position model (c./r.).
        """
        noncoding_point = self._noncoding.to_position(coordinate)

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

    def coordinate_to_coding(self, coordinate: int, degenerate: bool = False) -> CodingPoint:
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns CodingPoint: Coding position model (c./r.).
        """
        point = self._coordinate_to_coding(coordinate)

        region = point.region
        if not degenerate:
            return point

        if region == 'u':
            if self._coding[0] == 0:
                position = abs(point.offset)
            else:
                position = point.position + abs(point.offset)
            return CodingPoint(position=position, offset=0, region='-')
        if region == 'd':
            if self._exons[1] == self._coding[1]:
                position = abs(point.offset)
            else:
                position = point.position + abs(point.offset)
            return CodingPoint(position=position, offset=0, region='*')
        return point

    def coding_to_coordinate(self, point: CodingPoint) -> int:
        """Convert a coding position (c./r.) to a coordinate.

        :arg CodingPoint point: Coding position model (c./r.).

        :returns int: Coordinate.
        """

        region = point.region

        if region in ('u', 'd'):
            return self._noncoding.to_coordinate(point)

        position = point.position
        if region == '':
            noncoding_point = NonCodingPoint(
                position=position + self._coding[0] - 1,
                offset=point.offset,
                region=''
            )
            return self._noncoding.to_coordinate(noncoding_point)

        if region == '-':
            if position <= self._coding[0]:
                return self._noncoding.to_coordinate(
                    NonCodingPoint(
                        position=self._coding[0] - position,
                        offset=point.offset,
                        region=''
                    )
                )
            return self._noncoding.to_coordinate(
                NonCodingPoint(
                    position=0,
                    offset=point.offset + 1 - position,
                    region='u'
                )
            )

        return self._noncoding.to_coordinate(
            NonCodingPoint(
                position=self._coding[1] + position - 1,
                offset=point.offset,
                region=''
            )
        )

    def coordinate_to_protein(self, coordinate: int) -> ProteinPoint:
        """Convert a coordinate to a protein position (p.).

        :arg int coordinate: Coordinate.

        :returns ProteinPoint: Protein position model(p.).
        """
        point = self.coordinate_to_coding(coordinate)

        position = point.position
        if point.region in ('-', 'u'):
            return ProteinPoint(
                position=abs(-position // 3),
                position_in_codon=-position % 3 + 1,
                region=point.region,
                offset=point.offset
            )
        return ProteinPoint(
            position=(position + 2) // 3,
            position_in_codon=(position + 2) % 3 + 1,
            region=point.region,
            offset=point.offset
        )

    def protein_to_coordinate(self, point: ProteinPoint) -> int:
        """Convert a protein position (p.) to a coordinate.

        :arg ProteinPoint point: Protein position model(p.).

        :returns int: Coordinate.
        """
        if point.region in ('-', 'u'):
            return self.coding_to_coordinate(
                CodingPoint(
                    position=3 * point.position - point.position_in_codon + 1,
                    offset=point.offset,
                    region=point.region
                )
            )
        return self.coding_to_coordinate(
            CodingPoint(
                position=3 * point.position + point.position_in_codon - 3,
                offset=point.offset,
                region=point.region
            )
        )
