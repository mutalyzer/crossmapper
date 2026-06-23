from __future__ import annotations
from dataclasses import dataclass


def slotted_dataclass(cls=None, **kwargs):
    return dataclass(cls, slots=True, **kwargs)


# Basic dataclass module for locus and multi_locus
@slotted_dataclass
class Point:
    position: int
    offset: int = 0
    region: str = ''


@slotted_dataclass
class GenomicPoint:
    position: int

    def __post_init__(self) -> None:
        if not isinstance(self.position, int) or self.position <= 0:
            raise ValueError('Position must be a positive integer')

    def __str__(self) -> str:
        return f"{self.position}"


@slotted_dataclass
class NonCodingPoint(GenomicPoint):
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
            return f"{self.region}{self.position}"
        return f"{self.region}{self.position}{self.offset:+}"


@slotted_dataclass
class CodingPoint(NonCodingPoint):
    allowed_regions = {'', 'u', 'd', '-', '*'}


@slotted_dataclass
class ProteinPoint(CodingPoint):
    position_in_codon: int = 1

    def __post_init__(self) -> None:
        CodingPoint.__post_init__(self)

        if not isinstance(self.position_in_codon, int) or self.position_in_codon not in (1, 2, 3):
            raise ValueError('Position_in_codon must be 1, 2, or 3')
