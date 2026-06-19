from __future__ import annotations
from dataclasses import dataclass


# Basic dataclass module for locus and multi_locus
@dataclass(slots=True)
class Point:
    position: int
    offset: int = 0
    region: str = ''


@dataclass(slots=True)
class GenomicPoint:
    position: int

    def __post_init__(self) -> None:
        self._validate_position(self.position)

    def __str__(self) -> str:
        return f"{self.position}"

    @staticmethod
    def _validate_position(position: int) -> None:
        if not isinstance(position, int) or position < 0:
            raise ValueError('position must be a non-negative integer')


@dataclass(slots=True)
class NonCodingPoint(GenomicPoint):
    offset: int = 0
    region: str = ''

    allowed_regions = {'', 'u', 'd'}

    def __post_init__(self) -> None:
        GenomicPoint.__post_init__(self)
        self._validate_offset(self.offset)
        self._validate_region(self.region)

    @staticmethod
    def _validate_offset(offset: int) -> None:
        if not isinstance(offset, int):
            raise TypeError('offset must be an integer')

    def _validate_region(self, region: str) -> None:
        if not isinstance(region, str) or region not in self.allowed_regions:
            raise ValueError(f'region must be a string in {self.allowed_regions}')

    def __str__(self) -> str:
        if self.offset == 0:
            return f"{self.region}{self.position}"
        return f"{self.region}{self.position}{self.offset:+}"


@dataclass(slots=True)
class CodingPoint(NonCodingPoint):
    allowed_regions = {'', 'u', 'd', '-', '*'}


@dataclass(slots=True)
class ProteinPoint(CodingPoint):
    position_in_codon: int = 1

    def __post_init__(self) -> None:
        CodingPoint.__post_init__(self)
        self._validate_position_in_codon(self.position_in_codon)

    @staticmethod
    def _validate_position_in_codon(position_in_codon: int) -> None:
        if not isinstance(position_in_codon, int) or position_in_codon not in (1, 2, 3):
            raise ValueError('position_in_codon must be 1, 2, or 3')
