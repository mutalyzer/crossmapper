from mutalyzer_crossmapper.crossmapper import Coord, Genomic, NonCoding, Coding, GenomicPoint, NonCodingPoint, CodingPoint, ProteinPoint

from helper import degenerate_equal, invariant
import pytest

_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_cds = (32, 43)


def test_Genomic():
    """Genomic positions are coordinates incremented by one."""
    crossmap = Genomic()

    invariant(
        crossmap.coordinate_to_genomic,
        Coord(0),
        crossmap.genomic_to_coordinate,
        GenomicPoint(position=1),
    )
    invariant(
        crossmap.coordinate_to_genomic,
        Coord(98),
        crossmap.genomic_to_coordinate,
        GenomicPoint(position=99),
    )


def test_NonCoding():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons)

    # Boundary between upstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(3)    ,
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=1, offset=-2, region='u'),
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(4),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=1, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(5),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=1, offset=0, region=''),
    )

    # Boundary between downstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(71),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=22, offset=0, region=''),
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(72),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=22, offset=1, region='d'),
    )


def test_NonCoding_inverted():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons, inverted=True)

    # Boundary between upstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(72),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=1, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(71),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=1, offset=0, region=''),
    )

    # Boundary between downstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(5),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=22, offset=0, region=''),
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(4),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=22, offset=1, region='d'),
    )



def test_Coding():
    """Forward oriented coding transcript."""
    crossmap = Coding(_exons, _cds)

    # Boundary between 5' and CDS.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(31),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region='-'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(32),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region=''),
    )

    # Boundary between CDS and 3'.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(42),
        crossmap.coding_to_coordinate,
        CodingPoint(position=6, offset=0, region=''),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(43),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region='*'),
    )


def test_Coding_inverted():
    """Reverse oriented coding transcript."""
    crossmap = Coding(_exons, _cds, inverted=True)

    # Boundary between 5' and CDS.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(43),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region='-'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(42),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region=''),
    )

    # Boundary between CDS and 3'.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(32),
        crossmap.coding_to_coordinate,
        CodingPoint(position=6, offset=0, region=''),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(31),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region='*'),
    )


def test_Coding_regions():
    """The CDS can start or end on a region boundary."""
    crossmap = Coding([(10, 21), (30, 40), (49, 60)], (30, 40))

    # Upstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(25),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=5, region='-'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(26),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=-4, region=''),
    )

    # Downstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(44),
        crossmap.coding_to_coordinate,
        CodingPoint(position=10, offset=5, region=''),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(45),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=-4, region='*'),
    )


def test_Coding_regions_inverted():
    """The CDS can start or end on a region boundary."""
    crossmap = Coding([(10, 21), (30, 40), (49, 60)], (30, 40), inverted=True)

    # Upstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(44),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=5, region='-'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(43),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=-4, region=''),
    )

    # Downstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(25),
        crossmap.coding_to_coordinate,
        CodingPoint(position=10, offset=5, region=''),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(24),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=-4, region='*'),
    )


def test_Coding_no_utr5():
    """A 5' UTR may be missing."""
    crossmap = Coding([(10, 20)], (10, 15))

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(9),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(10),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region=''),
    )


def test_Coding_no_intron():
    crossmap = Coding([(10, 20), (20, 30)], (15, 25))

    invariant(
        crossmap.coordinate_to_coding,
        Coord(20),
        crossmap.coding_to_coordinate,
        CodingPoint(position=6, offset=0, region=''),
    )


def test_Coding_no_intron_inverted():
    crossmap = Coding([(10, 20), (20, 30)], (15, 25), inverted=True)

    invariant(
        crossmap.coordinate_to_coding,
        Coord(20),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=0, region=''),
    )


def test_Coding_one_base_intron():
    crossmap = Coding([(10, 19), (20, 30)], (15, 25))

    invariant(
        crossmap.coordinate_to_coding,
        Coord(19),
        crossmap.coding_to_coordinate,
        CodingPoint(position=4, offset=1, region=''),
    )


def test_Coding_one_base_intron_inverted():
    crossmap = Coding([(10, 19), (20, 30)], (15, 25), inverted=True)

    invariant(
        crossmap.coordinate_to_coding,
        Coord(19),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=1, region=''),
    )


def test_Coding_no_utr5_inverted():
    """A 5' UTR may be missing."""
    crossmap = Coding([(10, 20)], (15, 20), inverted=True)

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(20),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(19),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region=''),
    )


def test_Coding_no_utr3():
    """A 3' UTR may be missing."""
    crossmap = Coding([(10, 20)], (15, 20))
    print(crossmap.coordinate_to_coding(Coord(20)).position, crossmap.coordinate_to_coding(Coord(20)).offset, crossmap.coordinate_to_coding(Coord(20)).region)
    print(crossmap.coding_to_coordinate(CodingPoint(position=5, offset=1, region='d')))
    # Direct transition from CDS to downstream.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(19),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=0, region=''),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(20),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=1, region='d'),
    )


def test_Coding_no_utr3_inverted():
    """A 3' UTR may be missing."""
    crossmap = Coding([(10, 20)], (10, 15), inverted=True)

    # Direct transition from CDS to downstream.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(10),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=0, region=''),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(9),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=1, region='d'),
    )

def test_Coding_small_utr5():
    """A 5' UTR may be of length one."""
    crossmap = Coding([(10, 20)], (11, 15))

    # Transition from upstream to 5' UTR to CDS.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(9),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(10),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region='-'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(11),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region=''),
    )


def test_Coding_small_utr5_inverted():
    """A 5' UTR may be of length one."""
    crossmap = Coding([(10, 20)], (15, 19), inverted=True)

    # Transition from upstream to 5' UTR to CDS.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(20),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(19),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region='-'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(18),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region=''),
    )


def test_Coding_small_utr3():
    """A 5' UTR may be of length one."""
    crossmap = Coding([(10, 20)], (15, 19))

    # Transition from CDS to 3' UTR to downstream.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(18),
        crossmap.coding_to_coordinate,
        CodingPoint(position=4, offset=0, region=''),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(19),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region='*'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(20),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=1, region='d'),
    )


def test_Coding_small_utr3_inverted():
    """A 5' UTR may be of length one."""
    crossmap = Coding([(10, 20)], (11, 15), inverted=True)

    # Transition from CDS to 3' UTR to downstream.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(11),
        crossmap.coding_to_coordinate,
        CodingPoint(position=4, offset=0, region=''),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(10),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=0, region='*'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(9),
        crossmap.coding_to_coordinate,
        CodingPoint(position=1, offset=1, region='d'),
    )


def test_Coding_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    crossmap = Coding([(10, 20)], (11, 19))

    # Degenerate position in upstream.
    degenerate_equal(
        crossmap.coding_to_coordinate,
        Coord(9),
        [
            CodingPoint(position=1, offset=-1, region='u'),
            CodingPoint(position=2, offset=0, region='-'),
        ],
    )
    degenerate_equal(
        crossmap.coding_to_coordinate,
        Coord(20),
        [
            CodingPoint(position=1, offset=1, region='d'),
            CodingPoint(position=2, offset=0, region='*'),
        ],
    )


def test_Coding_inverted_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    crossmap = Coding([(10, 20)], (11, 19), inverted=True)

    degenerate_equal(
        crossmap.coding_to_coordinate,
        Coord(20),
        [
            CodingPoint(position=1, offset=-1, region='u'),
            CodingPoint(position=2, offset=0, region='-'),
        ],
    )
    degenerate_equal(
        crossmap.coding_to_coordinate,
        Coord(9),
        [
            CodingPoint(position=1, offset=1, region='d'),
            CodingPoint(position=2, offset=0, region='*'),
        ],
    )


def test_Coding_degenerate_return():
    """Degenerate upstream and downstream positions may be returned."""
    crossmap = Coding([(10, 20)], (11, 19))

    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=2, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(20), True) == CodingPoint(position=2, offset=0, region='*')


def test_Coding_inverted_degenerate_return():
    """Degenerate upstream and downstream positions may be returned."""
    crossmap = Coding([(10, 20)], (11, 19), inverted=True)

    assert crossmap.coordinate_to_coding(Coord(20), True) == CodingPoint(position=2, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(25), True) == CodingPoint(position=7, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=2, offset=0, region='*')


def test_Coding_no_utr5_degenerate_return():
    """A 5' UTR may be missing."""
    crossmap = Coding([(10, 20)], (10, 15))

    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=1, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(10), True) == CodingPoint(position=1, offset=0, region='')
    assert crossmap.coordinate_to_coding(Coord(19), True) == CodingPoint(position=5, offset=0, region='*')
    assert crossmap.coordinate_to_coding(Coord(20), True) == CodingPoint(position=6, offset=0, region='*')


def test_Coding_no_utr5_inverted_degenerate_return():
    """A 5' UTR may be missing."""
    crossmap = Coding([(10, 20)], (10, 15), inverted=True)

    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=1, offset=0, region='*')
    assert crossmap.coordinate_to_coding(Coord(10), True) == CodingPoint(position=5, offset=0, region='')
    assert crossmap.coordinate_to_coding(Coord(19), True) == CodingPoint(position=5, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(20), True) == CodingPoint(position=6, offset=0, region='-')


def test_Coding_no_utr3_degenerate_return():
    """A 3' UTR may be missing."""
    crossmap = Coding([(10, 20)], (15, 20))

    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=6, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(10), True) == CodingPoint(position=5, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(19), True) == CodingPoint(position=5, offset=0, region='')
    assert crossmap.coordinate_to_coding(Coord(20), True) == CodingPoint(position=1, offset=0, region='*')


def test_Coding_no_utr3_inverted_degenerate_return():
    """A 3' UTR may be missing."""

    crossmap = Coding([(10, 20)], (15, 20), inverted=True)

    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=6, offset=0, region='*')
    assert crossmap.coordinate_to_coding(Coord(10), True) == CodingPoint(position=5, offset=0, region='*')
    assert crossmap.coordinate_to_coding(Coord(19), True) == CodingPoint(position=1, offset=0, region='')
    assert crossmap.coordinate_to_coding(Coord(20), True) == CodingPoint(position=1, offset=0, region='-')


def test_Coding_small_utr5_degenerate_return():
    """A 5' UTR may be of length one."""
    crossmap = Coding([(10, 20)], (11, 15))

    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=2, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(10), True) == CodingPoint(position=1, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(11), True) == CodingPoint(position=1, offset=0, region='')


def test_Coding_small_utr5_inverted_degenerate_return():
    """A 5' UTR may be of length one."""
    crossmap = Coding([(10, 20)], (11, 15), inverted=True)

    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=2, offset=0, region='*')
    assert crossmap.coordinate_to_coding(Coord(10), True) == CodingPoint(position=1, offset=0, region='*')
    assert crossmap.coordinate_to_coding(Coord(11), True) == CodingPoint(position=4, offset=0, region='')


def test_Coding_small_utr3_degenerate_return():
    """A 3' UTR may be of length one."""
    crossmap = Coding([(10, 20)], (15, 19), inverted=False)

    assert crossmap.coordinate_to_coding(Coord(18), True) == CodingPoint(position=4, offset=0, region='')
    assert crossmap.coordinate_to_coding(Coord(19), True) == CodingPoint(position=1, offset=0, region='*')
    assert crossmap.coordinate_to_coding(Coord(20), True) == CodingPoint(position=2, offset=0, region='*')


def test_Coding_small_utr3_inverted_degenerate_return():
    """A 3' UTR may be of length one."""
    crossmap = Coding([(10, 20)], (15, 19), inverted=True)

    assert crossmap.coordinate_to_coding(Coord(18), True) == CodingPoint(position=1, offset=0, region='')
    assert crossmap.coordinate_to_coding(Coord(19), True) == CodingPoint(position=1, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(20), True) == CodingPoint(position=2, offset=0, region='-')


def test_Coding_two_exons_inverted_degenerate_return():
    """Degenerate upstream and downstream positions may be returned."""
    crossmap = Coding([(10, 20), (30, 40)], (18, 37), inverted=True)

    assert crossmap.coordinate_to_coding(Coord(5), True) == CodingPoint(position=13, offset=0, region='*')
    assert crossmap.coordinate_to_coding(Coord(25), True) == CodingPoint(position=7, offset=5, region='')
    assert crossmap.coordinate_to_coding(Coord(35), True) == CodingPoint(position=2, offset=0, region='')
    assert crossmap.coordinate_to_coding(Coord(38), True) == CodingPoint(position=2, offset=0, region='-')


def test_Coding_degenerate_no_return():
    """Degenerate internal positions do not exist."""
    crossmap = Coding([(10, 20), (30, 40)], (10, 40))

    assert crossmap.coordinate_to_coding(Coord(25)) == crossmap.coordinate_to_coding(Coord(25), True)


def test_Coding_inverted_degenerate_no_return():
    """Degenerate internal positions do not exist."""
    crossmap = Coding([(10, 20), (30, 40)], (10, 40), inverted=True)

    assert crossmap.coordinate_to_coding(Coord(25)) == crossmap.coordinate_to_coding(Coord(25), True)


def test_Coding_no_utr_degenerate():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11))

    degenerate_equal(
        crossmap.coding_to_coordinate,
        Coord(9),
        [
            CodingPoint(position=1, offset=-1, region='u'),
            CodingPoint(position=1, offset=0, region='-'),
        ],
    )
    degenerate_equal(
        crossmap.coding_to_coordinate,
        Coord(11),
        [
            CodingPoint(position=1, offset=1, region='d'),
            CodingPoint(position=1, offset=0, region='*'),
        ],
    )


def test_Coding_inverted_no_utr_degenerate():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11), inverted=True)

    degenerate_equal(
        crossmap.coding_to_coordinate,
        Coord(11),
        [
            CodingPoint(position=1, offset=-1, region='u'),
            CodingPoint(position=1, offset=0, region='-'),
            # CodingPoint(position=1, offset=-2, region='*'),
            # CodingPoint(position=1, offset=-1, region=''),
            # CodingPoint(position=1, offset=-1, region='d'),
        ],
    )
    degenerate_equal(
        crossmap.coding_to_coordinate,
        Coord(9),
        [
            CodingPoint(position=1, offset=1, region='d'),
            CodingPoint(position=1, offset=0, region='*'),
            # CodingPoint(position=1, offset=2, region='-'),
            # CodingPoint(position=1, offset=1, region=''),
            # CodingPoint(position=1, offset=1, region='u'),
        ],
    )


def test_Coding_no_utr_degenerate_return():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11))

    assert crossmap.coordinate_to_coding(Coord(8), True) == CodingPoint(position=2, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=1, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(11), True) == CodingPoint(position=1, offset=0, region='*')
    assert crossmap.coordinate_to_coding(Coord(12), True) == CodingPoint(position=2, offset=0, region='*')


def test_Coding_inverted_no_utr_degenerate_return():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11), inverted=True)

    assert crossmap.coordinate_to_coding(Coord(11), True) == CodingPoint(position=1, offset=0, region='-')
    assert crossmap.coordinate_to_coding(Coord(9), True) == CodingPoint(position=1, offset=0, region='*')


def test_Coding_protein():
    """Protein positions."""
    crossmap = Coding(_exons, _cds)
    for i in range(0, 80):
        t = crossmap.coordinate_to_protein(Coord(i))
        print(f"{i}: {crossmap.coordinate_to_coding(Coord(i))},{t.region} {t.position} {t.offset} {t.position_in_codon}")

    print(crossmap.coordinate_to_protein(Coord(4)))

    print(crossmap.protein_to_coordinate(ProteinPoint(position=4, offset=-1, region='u', position_in_codon=2)))

    # Boundary between upstream and 5' UTR
    invariant(
        crossmap.coordinate_to_protein,
        Coord(4),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=4, offset=-1, region='u', position_in_codon=2)
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(5),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=4, offset=0, region='-', position_in_codon=2)
    )

    # Boundary between 5' UTR and CDS
    invariant(
        crossmap.coordinate_to_protein,
        Coord(31),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=1, offset=0, region='-', position_in_codon=3),
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(32),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=1, offset=0, region='', position_in_codon=1),
    )

    # Intron boundary.
    invariant(
        crossmap.coordinate_to_protein,
        Coord(34),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=1, offset=0, region='', position_in_codon=3),
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(35),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=1, offset=1, region='', position_in_codon=3),
    )

    # Boundary between CDS and 3' UTR.
    invariant(
        crossmap.coordinate_to_protein,
        Coord(42),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=2, offset=0, region='', position_in_codon=3),
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(43),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=1, offset=0, region='*', position_in_codon=1),
    )

    # Boundary between 3' UTR and downstream
    invariant(
        crossmap.coordinate_to_protein,
        Coord(71),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=2, offset=0, region='*', position_in_codon=2)
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(72),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=2, offset=1, region='d', position_in_codon=2)
    )


def test_Coding_inverted_protein():
    """Protein positions."""
    crossmap = Coding(_exons, _cds, inverted = True)

    # Boundary between upstream and 5' UTR
    invariant(
        crossmap.coordinate_to_protein,
        Coord(4),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=4, offset=1, region='d', position_in_codon=2)
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(5),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=4, offset=0, region='*', position_in_codon=2)
    )

    # Boundary between 5' UTR and CDS
    invariant(
        crossmap.coordinate_to_protein,
        Coord(31),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=1, offset=0, region='*', position_in_codon=1),
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(32),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=2, offset=0, region='', position_in_codon=3),
    )

    # Intron boundary.
    invariant(
        crossmap.coordinate_to_protein,
        Coord(34),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=2, offset=0, region='', position_in_codon=1),
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(35),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=2, offset=-1, region='', position_in_codon=1),
    )

    # Boundary between CDS and 3' UTR.
    invariant(
        crossmap.coordinate_to_protein,
        Coord(42),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=1, offset=0, region='', position_in_codon=1),
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(43),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=1, offset=0, region='-', position_in_codon=3),
    )

    # Boundary between 3' UTR and downstream
    invariant(
        crossmap.coordinate_to_protein,
        Coord(71),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=2, offset=0, region='-', position_in_codon=2)
    )
    invariant(
        crossmap.coordinate_to_protein,
        Coord(72),
        crossmap.protein_to_coordinate,
        ProteinPoint(position=2, offset=-1, region='u', position_in_codon=2)
    )
