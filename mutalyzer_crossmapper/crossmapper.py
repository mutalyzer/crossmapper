from dataclasses import dataclass

from .multi_locus import MultiLocus, Point, _check_in_range, _check_multi_locus
from .locus import Coord, _check_locus, _check_int
from .location import _nearest_location

@dataclass(slots=True)
class GenomicPoint:
    """Genomic dataclass."""
    position: int

    def __post_init__(self) -> None:
        _check_int(self.position)
        if self.position <= 0:
            raise ValueError(f'Position {self.position} must be a positive integer.')

    def __str__(self) -> str:
        return f'{self.position}'


class Genomic():
    """Genomic crossmap object."""
    def coordinate_to_genomic(self, coord: Coord, length: int | None = None) -> GenomicPoint:
        """Convert a coordinate dataclass to a genomic point dataclass (g./m./o.).

        :arg Coord coordinate: Coordinate dataclass.
        :arg int|None length: Length of the sequence.

        :returns GenomicPoint: Genomic point dataclass.
        """
        if length is not None:
            _check_in_range(coord.coordinate, length)
        return GenomicPoint(coord.coordinate + 1)

    def genomic_to_coordinate(self, point: GenomicPoint) -> Coord:
        """Convert a genomic point dataclass (g./m./o.) to a coordinate dataclass.

        :arg GenomicPoint point: Genomic point dataclass.

        :returns Coord: Coordinate dataclass.
        """
        return Coord(point.position - 1)


@dataclass(slots=True)
class NonCodingPoint(GenomicPoint):
    """NonCoding dataclass."""
    offset: int = 0
    region: str = ''

    allowed_regions = ['', 'u', 'd']

    def __post_init__(self) -> None:
        # Python version 3.11 and 3.10: cannot use super() due to conflicts with slots=True
        GenomicPoint.__post_init__(self)

        _check_int(self.offset)
        if self.region not in self.allowed_regions:
            raise ValueError(
                f'Region {self.region} is invalid, it must be a string from {self.allowed_regions}.'
            )

    def __str__(self) -> str:
        if self.offset == 0:
            return f'{self.region}{self.position}'
        if self.region in ('u', 'd'):
            return f'{self.region}{abs(self.offset)}'
        return f'{self.region}{self.position}{self.offset:+}'


class NonCoding(Genomic):
    """NonCoding crossmap object."""

    def __init__(
        self,
        locations: list[tuple[int, int]],
        length: int | None = None,
        inverted: bool = False,
    ) -> None:
        """
        :arg list locations: List of locus locations.
        :arg int|None length: Length of the reference sequence.
        :arg bool inverted: Orientation.
        """
        _check_multi_locus(locations, length)
        self._inverted = inverted
        self._noncoding = MultiLocus(locations, length, inverted)

    def coordinate_to_noncoding(self, coord: Coord) -> NonCodingPoint:
        """Convert a coordinate dataclass to a noncoding point dataclass (n./r.).

        :arg Coord coord: Coordinate dataclass.

        :returns NonCodingPoint: Noncoding point dataclass.
        """
        point = self._noncoding.to_position(coord)
        return NonCodingPoint(
            position=point.position + 1,
            offset=point.offset,
            region=point.region
        )

    def noncoding_to_coordinate(self, point: NonCodingPoint) -> Coord:
        """Convert a noncoding point dataclass (n./r.) to a coordinate dataclass.

        :arg NonCodingPoint point: Noncoding point dataclass.

        :returns Coord: Coordinate dataclass.
        """
        try:
            return self._noncoding.to_coordinate(
                Point(
                    position=point.position - 1,
                    offset=point.offset,
                    region=point.region
                )
            )
        except ValueError as error:
            if 'Position' in str(error):
                raise ValueError(str(error).replace(str(point.position - 1), str(point.position)))
            raise error
        except IndexError as error:
            if 'Position' in str(error):
                raise IndexError(str(error).replace(str(point.position - 1), str(point.position)))
            raise error


@dataclass(slots=True)
class CodingPoint(NonCodingPoint):
    """Coding dataclass."""
    allowed_regions = ['', 'u', 'd', '-', '*']


@dataclass(slots=True)
class ProteinPoint(CodingPoint):
    """Protein dataclass."""
    position_in_codon: int = 1

    def __post_init__(self) -> None:
        CodingPoint.__post_init__(self)

        if not isinstance(self.position_in_codon, int) or self.position_in_codon not in (1, 2, 3):
            raise ValueError('Position_in_codon must be 1, 2, or 3.')

    def __str__(self) -> str:
        if self.offset == 0 and self.region == '':
            return f'{self.position}'
        return '?'


class Coding(NonCoding):
    """Coding crossmap object."""
    def __init__(
            self,
            locations: list[tuple[int, int]],
            cds: tuple[int, int],
            length: int|None = None,
            inverted: bool = False
    ) -> None:
        """
        :arg list locations: List of locus locations.
        :arg tuple cds: Locus location.
        :arg int|None length: Length of the reference sequence.
        :arg bool inverted: Orientation.
        """
        NonCoding.__init__(self, locations, length, inverted)
        self._check_cds(cds, locations, length)

        cds_start = self._noncoding.to_position(Coord(cds[0]))
        cds_end = self._noncoding.to_position(Coord(cds[1] - 1))
        exon_start = self._noncoding.to_position(Coord(locations[0][0]))
        exon_end = self._noncoding.to_position(Coord(locations[-1][1] - 1))

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

    def _check_cds(
        self,
        cds: tuple[int, int],
        locations: list[tuple[int, int]],
        length: int | None = None,
    ) -> None:
        """Check if the CDS is valid."""
        _check_locus(cds)
        if length is not None:
            _check_in_range(cds[1], length)
        for coord in cds:
            index = _nearest_location(locations, coord)
            if coord < locations[index][0] or coord > locations[index][1]:
                raise ValueError(f'Coordinate {coord} of CDS {cds} is not within any exon.')

    def _validate_point(self, position: int, region: str) -> None:
        """Validate a coding point model under HGVS rules.

        :arg int position: Position.
        :arg str region: Region.
        """
        if region == 'u':
            if position not in (1, self._coding[0]):
                raise ValueError(f'Position {position} is not in upstream boundary.')
        if region == '-':
            if position not in range(1, self._coding[0] + 1):
                raise ValueError(f'Position {position} exceeds - region.')
        if region == '':
            if position not in range(1, self._coding[1] - self._coding[0] + 1):
                raise ValueError(f'Position {position} exceeds coding region.')
        if region == '*':
            if position not in range(1, self._exons[1] - self._coding[1] + 1):
                raise ValueError(f'Position {position} exceeds * region.')
        if region == 'd':
            if position not in (1, self._coding[0], self._exons[1] - self._coding[1]):
                raise ValueError(f'Position {position} is not in downstream boundary.')


    def _coordinate_to_coding(self, coord: Coord) -> CodingPoint:
        """Convert a coordinate dataclass to a coding point dataclass (c./r.).

        :arg Coord coord: Coordinate dataclass.

        :returns CodingPoint: Coding point dataclass (c./r.).
        """
        noncoding_point = self._noncoding.to_position(coord)

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

    def coordinate_to_coding(self, coord: Coord, degenerate: bool = False) -> CodingPoint:
        """Convert a coordinate dataclass to a coding point dataclass (c./r.).

        :arg Coord coord: Coordinate dataclass.
        :arg bool degenerate: Return a degenerate coding point dataclass.

        :returns CodingPoint: Coding point dataclass (c./r.).
        """
        point = self._coordinate_to_coding(coord)

        if not degenerate:
            return point

        offset = point.offset
        if point.region == 'u':
            if self._coding[0] == 0:
                position = -offset
            else:
                position = point.position - offset
            return CodingPoint(position=position, offset=0, region='-')
        if point.region == 'd':
            if self._exons[1] == self._coding[1]:
                position = offset
            else:
                position = point.position + offset
            return CodingPoint(position=position, offset=0, region='*')
        return point

    def _coding_to_coordinate(self, point: CodingPoint) -> Coord:
        """Convert a coding point dataclass (c./r.) to a coordinate dataclass.

        :arg CodingPoint point: Coding point dataclass (c./r.).

        :returns Coord: Coordinate dataclass.
        """
        region = point.region
        position = point.position

        self._validate_point(position, region)

        # For missing 3' UTR or 5' UTR
        if region in ('u', 'd'):
            if region == 'u':
                position = 1
            if region == 'd':
                if self._coding[1] == self._exons[1]:
                    position = self._coding[1]
                else:
                    position = position + self._coding[1]
            return self._noncoding.to_coordinate(
                Point(position=position - 1, region=point.region, offset=point.offset)
            )

        if region == '':
            position = position + self._coding[0] - 1
        elif region == '-':
            position = self._coding[0] - position
        elif region == '*':
            position = self._coding[1] + position - 1

        try:
            return self._noncoding.to_coordinate(
                Point(position=position, region='', offset=point.offset)
            )
        except ValueError as error:
            if 'Position' in str(error):
                raise ValueError(str(error).replace(str(position), str(point.position)))
            raise

    def coding_to_coordinate(self, point: CodingPoint) -> Coord:
        """Convert a coding point dataclass (c./r.) to a coordinate dataclass.

        :arg CodingPoint point: Coding point dataclass (c./r.).

        :returns Coord: Coordinate dataclass.
        """
        # Silently correct for degenerate points
        if point.offset == 0:
            if point.region == '-' and point.position > self._coding[0]:
                if self._coding[0] == 0:
                    return self._coding_to_coordinate(CodingPoint(
                        position=1,
                        offset=self._coding[0] - point.position,
                        region='u',
                    ))
                return self._coding_to_coordinate(CodingPoint(
                    position=self._coding[0],
                    offset=self._coding[0] - point.position,
                    region='u',
                ))
            if point.region == '*' and point.position > self._exons[1] - self._coding[1]:
                if self._exons[1] == self._coding[1]:
                    return self._coding_to_coordinate(CodingPoint(
                        position=1,
                        offset=point.position - (self._exons[1] - self._coding[1]),
                        region='d',
                    ))
                return self._coding_to_coordinate(CodingPoint(
                    position=self._exons[1] - self._coding[1],
                    offset=point.position - (self._exons[1] - self._coding[1]),
                    region='d',
                ))

        return self._coding_to_coordinate(point)

    def coordinate_to_protein(self, coord: Coord) -> ProteinPoint:
        """Convert a coordinate dataclass to a protein point dataclass (p.).

        :arg Coord coord: Coordinate dataclass.

        :returns ProteinPoint: Protein point dataclass (p.).
        """
        point = self.coordinate_to_coding(coord)

        position = point.position
        if point.region in ('-', 'u'):
            return ProteinPoint(
                position=-(-position // 3),
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

    def protein_to_coordinate(self, point: ProteinPoint) -> Coord:
        """Convert a protein point dataclass (p.) to a coordinate dataclass.

        :arg ProteinPoint point: Protein point dataclass (p.).

        :returns Coord: Coordinate dataclass.
        """
        if point.region in ('-', 'u'):
            return self._coding_to_coordinate(
                CodingPoint(
                    position=3 * point.position - point.position_in_codon + 1,
                    offset=point.offset,
                    region=point.region
                )
            )
        return self._coding_to_coordinate(
            CodingPoint(
                position=3 * point.position + point.position_in_codon - 3,
                offset=point.offset,
                region=point.region
            )
        )
