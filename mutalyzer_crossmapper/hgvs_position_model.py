"""
HGVS Position Model -
    a dataclass object to bridge HGVS position component and Crossmapper outputs.
"""
from dataclasses import dataclass
from typing import Optional

@dataclass
class HGVSPositionModel:
    """
    Represent the position component of an HGVS variant description.
    This model captures details necessary to describe the '[position]' part in an HGVS
    description of the form
        [reference sequence]:[sequence type].[position][variant type][change]
    """
    position: int
    offset: Optional[int] = None
    region: Optional[str] = None
    position_in_codon: Optional[int] = None


    def __post_init__(self):
        # validate position
        if self.position <= 0:
            raise ValueError("Position must be a positive integer.")

        # validate region
        region_values = {"u", "-", "", "*", "d"}
        if self.region is not None and self.region not in region_values:
            raise ValueError(
                f"Invalid region value: {self.region}. Allowed values are: {region_values}"
            )

        # validate position_in_codon
        codon_values = {1, 2, 3}
        if self.position_in_codon is not None and self.position_in_codon not in codon_values:
            raise ValueError(
                f"Invalid position in codon value: {self.position_in_codon}. "
                f"Allowed values are: {codon_values}"
            )
