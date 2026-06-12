from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

NONCODING_REGIONS = {'', 'u', 'd'}
CODING_REGIONS = NONCODING_REGIONS | {'-', '*'}


@dataclass(slots=True)
class GenomicPoint:
    position: int

    def __post_init__(self) -> None:
        self._validate_position(self.position)

    @staticmethod
    def _validate_position(position: int) -> None:
        if not isinstance(position, int) or position < 0:
            raise TypeError("position must be a non-negative integer")

    def to_dict(self) -> dict[str, Any]:
        return {'position': self.position}

    @classmethod
    def to_dataclass(cls, point: Any) -> GenomicPoint:
        if isinstance(point, cls):
            return point
        if isinstance(point, dict):
            return cls(position=point["position"])
        raise TypeError(f"Cannot convert {type(point)}")


@dataclass(slots=True)
class NonCodingPoint(GenomicPoint):
    offset: int = 0
    region: str = ""

    allowed_regions = NONCODING_REGIONS

    def __post_init__(self) -> None:
        GenomicPoint.__post_init__(self)
        self._validate_offset(self.offset)
        self._validate_region(self.region)

    @staticmethod
    def _validate_offset(offset: int) -> None:
        if not isinstance(offset, int):
            raise TypeError("offset must be an integer")

    def _validate_region(self, region: str) -> None:
        if not isinstance(region, str) or region not in self.allowed_regions:
            raise ValueError(f"region must be a string in {self.allowed_regions}")

    def to_dict(self) -> dict[str, Any]:
        return {
            'position': self.position,
            'offset': self.offset,
            'region': self.region,
        }

    @classmethod
    def to_dataclass(cls, point: Any) -> NonCodingPoint:
        if isinstance(point, cls):
            return point

        if isinstance(point, GenomicPoint):
            return cls(position=point.position)

        if isinstance(point, dict):
            return cls(
                position=point["position"],
                offset=point.get("offset", 0),
                region=point.get("region", ""),
            )

        raise TypeError(f"Cannot convert {type(point)}")


@dataclass(slots=True)
class CodingPoint(NonCodingPoint):
    allowed_regions = CODING_REGIONS

    @classmethod
    def to_dataclass(cls, point: Any) -> CodingPoint:
        if isinstance(point, cls):
            return point

        if isinstance(point, NonCodingPoint):
            return cls(
                position=point.position,
                offset=point.offset,
                region=point.region,
            )

        if isinstance(point, dict):
            return cls(
                position=point["position"],
                offset=point.get("offset", 0),
                region=point.get("region", ""),
            )

        raise TypeError(f"Cannot convert {type(point)}")


@dataclass(slots=True)
class ProteinPoint(CodingPoint):
    position_in_codon: int = 1

    def __post_init__(self) -> None:
        CodingPoint.__post_init__(self)
        self._validate_position_in_codon(self.position_in_codon)

    def to_dict(self) -> dict[str, Any]:
        return {
            'position': self.position,
            'offset': self.offset,
            'region': self.region,
            'position_in_codon': self.position_in_codon,
        }

    @staticmethod
    def _validate_position_in_codon(position_in_codon: int) -> None:
        if not isinstance(position_in_codon, int) or position_in_codon not in (1, 2, 3):
            raise ValueError("position_in_codon must be 1, 2, or 3")

    @classmethod
    def to_dataclass(cls, point: Any) -> ProteinPoint:
        if isinstance(point, cls):
            return point

        if isinstance(point, CodingPoint):
            return cls(
                position=point.position,
                offset=point.offset,
                region=point.region,
                position_in_codon=getattr(point, "position_in_codon", 1),
            )

        if isinstance(point, dict):
            return cls(
                position=point["position"],
                offset=point.get("offset", 0),
                region=point.get("region", ""),
                position_in_codon=point.get("position_in_codon", 1),
            )

        raise TypeError(f"Cannot convert {type(point)}")


