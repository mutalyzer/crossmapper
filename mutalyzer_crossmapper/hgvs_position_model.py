"""
HGVS Position Model -
    a dataclass object to bridge HGVS position component and Crossmapper outputs.
"""
from dataclasses import dataclass
from typing import Optional, Tuple

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


    # Convert from tuple to HGVSPositionModel

    #TODO: check for inverted and degerate options, now only support non-inverted and non-degenerate cases
    @classmethod
    def to_hgvs_position_model(cls, raw_tuple:Tuple):
        """Convert crossmapper tuple to an HGVSPositionModel instance."""
        if not raw_tuple:
            raise ValueError("Input tuple position cannot be empty.")

        # Genomic
        if len(raw_tuple) == 1:
            return cls(position=raw_tuple[0])
        # Non-coding
        if len(raw_tuple) == 3:
            pass

        # Coding
        #(c_pos, offset, in_cds, offset_to_exon_boundary)
        if len(raw_tuple) == 4:
            c_pos, offset, cds, dis_to_exon_boundary = raw_tuple
            region = cls._determine_region(cds, dis_to_exon_boundary)
            return cls(position=c_pos, offset=offset, region=region)

        # Protein (
        if len(raw_tuple) == 5:
            p_pos, codon_pos, offset, cds, dis_to_exon_boundary = raw_tuple
            if cds == 0: # in CDS
                return cls(
                    position=p_pos,
                    region="",
                    position_in_codon=codon_pos
                )
            else:
                # TODO: shall we support HGVSPositionModel outside of CDS for protein?
                pass


    @staticmethod
    def _determine_region(cds, dis_to_exon_boundary):
        if dis_to_exon_boundary < 0:
            return "u"
        elif dis_to_exon_boundary > 0:
            return "d"
        else: # in translation range, check if in CDS or not
            if cds < 0: # before CDS
                return "-"
            elif cds > 0: # after CDS
                return "*"
            else:
                return ""