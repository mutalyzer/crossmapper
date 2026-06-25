from dataclasses import dataclass

from .multi_locus import MultiLocus
from .locus import Point

@dataclass(slots=True)
class GenomicPoint:
    """Genomic dataclass."""
    position: int

    def __post_init__(self) -> None:
        if not isinstance(self.position, int) or self.position <= 0:
            raise ValueError('Position must be a positive integer')

    def __str__(self) -> str:
        return f'{self.position}'


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


@dataclass(slots=True)
class NonCodingPoint(GenomicPoint):
    """NonCoding dataclass."""
    offset: int = 0
    region: str = ''

    allowed_regions = {'', 'u', 'd'}

    def __post_init__(self) -> None:
        # Python version 3.11 and 3.10: cannot use super() due to conflicts with slots=True
        GenomicPoint.__post_init__(self)

        if not isinstance(self.offset, int):
            raise TypeError('Offset must be an integer')
        if not isinstance(self.region, str) or self.region not in self.allowed_regions:
            raise ValueError(f'Region must be a string in {self.allowed_regions}')

    def __str__(self) -> str:
        if self.offset == 0:
            return f'{self.region}{self.position}'
        return f'{self.region}{self.position}{self.offset:+}'


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
            Point(
                position=point.position - 1,
                offset=point.offset,
                region=point.region
            )
        )


@dataclass(slots=True)
class CodingPoint(NonCodingPoint):
    """Coding dataclass."""
    allowed_regions = {'', 'u', 'd', '-', '*'}


@dataclass(slots=True)
class ProteinPoint(CodingPoint):
    """Protein dataclass."""
    position_in_codon: int = 1

    def __post_init__(self) -> None:
        CodingPoint.__post_init__(self)

        if not isinstance(self.position_in_codon, int) or self.position_in_codon not in (1, 2, 3):
            raise ValueError('Position_in_codon must be 1, 2, or 3')


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
        """Convert a coordinate to a coding point model (c./r.).

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
        """Convert a coordinate to a coding point model (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns CodingPoint: Coding point model (c./r.).
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

        :arg CodingPoint point: Coding point model (c./r.).

        :returns int: Coordinate.
        """
        region = point.region

        if region in ('u', 'd'):
            return self._noncoding.to_coordinate(
                Point(position=point.position, region=point.region, offset=point.offset)
            )

        position = point.position
        if region == '':
            position = position + self._coding[0] - 1
        elif region == '-':
            position = self._coding[0] - position
        elif region == '*':
            position = self._coding[1] + position - 1
        return self._noncoding.to_coordinate(
            Point(position=position, region='', offset=point.offset)
        )

    def coordinate_to_protein(self, coordinate: int) -> ProteinPoint:
        """Convert a coordinate to a protein point model (p.).

        :arg int coordinate: Coordinate.

        :returns ProteinPoint: Protein point model(p.).
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

        :arg ProteinPoint point: Protein point model(p.).

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
