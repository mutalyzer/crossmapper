from mutalyzer_crossmapper.crossmapper import Coord, Genomic, NonCoding, Coding, GenomicPoint, NonCodingPoint, CodingPoint, ProteinPoint

from helper import degenerate_equal, invariant
import pytest

_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_cds = (32, 43)


def test_GenomicPoint_invalid_initialization():
    """GenomicPoint cannot be initialized with invalid position."""
    with pytest.raises(ValueError) as e:
        GenomicPoint(position=0)
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        GenomicPoint(position=-1)
    assert str(e.value) == "Position -1 must be a positive integer."
    with pytest.raises(ValueError) as e:
        GenomicPoint(position=[101])
    assert str(e.value) == "Value must be an integer."


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


def test_Genomic_invalid_with_length():
    """Raise ValueError if coordinate is out of bounds."""
    crossmap = Genomic()
    with pytest.raises(ValueError) as e:
        crossmap.coordinate_to_genomic(Coord(-1), 99)
    assert str(e.value) == "Value must be non-negative."
    with pytest.raises(ValueError) as e:
        crossmap.coordinate_to_genomic(Coord(99), 99)
    assert str(e.value) == "Value 99 must be within the bounds of the reference length 99."


def test_Genomic_with_length():
    """Genomic positions are coordinates incremented by one."""
    crossmap = Genomic()

    invariant(
        crossmap.coordinate_to_genomic,
        (Coord(0), 99),
        crossmap.genomic_to_coordinate,
        GenomicPoint(position=1),
    )
    invariant(
        crossmap.coordinate_to_genomic,
        (Coord(98), 99),
        crossmap.genomic_to_coordinate,
        GenomicPoint(position=99),
    )


def test_NonCodingPoint_invalid_initialization():
    """Raise error with invalid initialization."""
    with pytest.raises(ValueError) as e:
        NonCodingPoint(position=0, offset=0, region='u')
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        NonCodingPoint(position=0, offset=0, region='d')
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        NonCodingPoint(position=0, offset=0, region='')
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        NonCodingPoint(position=0, offset=0, region='*')
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        NonCodingPoint(position=-1, offset=0, region='')
    assert str(e.value) == "Position -1 must be a positive integer."
    with pytest.raises(ValueError) as e:
        NonCodingPoint(position=1, offset=None, region='u')
    assert str(e.value) == "Value must be an integer."
    with pytest.raises(ValueError) as e:
        NonCodingPoint(position=1, offset=1, region='-')
    assert str(e.value) == "Region must be a string in ['', 'u', 'd']."


def test_NonCoding_invalid():
    """Raise ValueError if noncoding is invalid."""
    with pytest.raises(ValueError) as e:
        NonCoding([()])
    assert str(e.value) == "Locus must be a tuple of two values."
    with pytest.raises(ValueError) as e:
        NonCoding([(10)])
    assert str(e.value) == "Locus must be a tuple of two values."
    with pytest.raises(ValueError) as e:
        NonCoding([(10, 20), (15, 25)])
    assert str(e.value) == "Locus (15, 25) and locus (10, 20) are overlapping."
    with pytest.raises(ValueError) as e:
        NonCoding([(None, 20), (30, None)])
    assert str(e.value) == "Value must be an integer."
    with pytest.raises(ValueError) as e:
        NonCoding(_exons, length=70)
    assert str(e.value) == "Value 72 must be within the bounds of the reference length 70."

    # Reverse orientation
    with pytest.raises(ValueError) as e:
        NonCoding([()], inverted=True)
    assert str(e.value) == "Locus must be a tuple of two values."
    with pytest.raises(ValueError) as e:
        NonCoding([(10)], inverted=True)
    assert str(e.value) == "Locus must be a tuple of two values."
    with pytest.raises(ValueError) as e:
        NonCoding([(10, 20), (15, 25)], inverted=True)
    assert str(e.value) == "Locus (15, 25) and locus (10, 20) are overlapping."
    with pytest.raises(ValueError) as e:
        NonCoding([(None, 20), (30, None)], inverted=True)
    assert str(e.value) == "Value must be an integer."
    with pytest.raises(ValueError) as e:
        NonCoding(_exons, length=70, inverted=True)
    assert str(e.value) == "Value 72 must be within the bounds of the reference length 70."


def test_NonCoding_invalid_with_length():
    """Raise ValueError if coordinate is out of bounds."""
    with pytest.raises(ValueError) as e:
        NonCoding(_exons, length=70)
    assert str(e.value) == "Value 72 must be within the bounds of the reference length 70."
    # Reverse orientation
    with pytest.raises(ValueError) as e:
        NonCoding(_exons, length=70, inverted=True)
    assert str(e.value) == "Value 72 must be within the bounds of the reference length 70."


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


def test_NonCoding_with_length():
    """Raise ValueError if coordinate is out of bounds."""
    crossmap = NonCoding(_exons, length=75)

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

    # Boundary between downstream and sequence end.
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(74),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=22, offset=3, region='d'),
    )
    with pytest.raises(ValueError) as e:
        crossmap.coordinate_to_noncoding(Coord(75))
    assert str(e.value) == "Value 75 must be within the bounds of the reference length 75."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=22, offset=4, region='d'))
    assert str(e.value) == "Offset 4 exceeds downstream region."


def test_NonCoding_inverted():
    """Reverse oriented noncoding transcript."""
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


def test_NonCoding_inverted_with_length():
    """Reverse oriented noncoding transcript."""
    crossmap = NonCoding(_exons, length=75, inverted=True)

    # Boundary between upstream and sequence end.
    with pytest.raises(ValueError) as e:
        crossmap.coordinate_to_noncoding(Coord(75))
    assert str(e.value) == "Value 75 must be within the bounds of the reference length 75."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=-4, region='u'))
    assert str(e.value) == "Offset -4 exceeds upstream region."
    invariant(
        crossmap.coordinate_to_noncoding,
        Coord(74),
        crossmap.noncoding_to_coordinate,
        NonCodingPoint(position=1, offset=-3, region='u'),
    )

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


def test_NonCoding_invalid_position():
    """Raise error if position is not valid under HGVS rules."""
    crossmap = NonCoding(_exons, length=75)
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=0, offset=1, region='u'))
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=2, offset=1, region='u'))
    assert str(e.value) == "Position 2 is not at upstream boundary."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=23, offset=0, region=''))
    assert str(e.value) == "Position 23 exceeds multi locus length."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=21, offset=1, region='d'))
    assert str(e.value) == "Position 21 is not at downstream boundary."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=30, offset=-1, region='d'))
    assert str(e.value) == "Position 30 is not at downstream boundary."


def test_NonCoding_invalid_position_inverted():
    """Raise error if position is not valid under HGVS rules."""
    crossmap = NonCoding(_exons, length=75, inverted=True)
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=0, offset=1, region='u'))
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=2, offset=1, region='u'))
    assert str(e.value) == "Position 2 is not at upstream boundary."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=23, offset=0, region=''))
    assert str(e.value) == "Position 23 exceeds multi locus length."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=21, offset=1, region='d'))
    assert str(e.value) == "Position 21 is not at downstream boundary."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=30, offset=-1, region='d'))
    assert str(e.value) == "Position 30 is not at downstream boundary."


def test_NonCoding_invalid_offset():
    """Raise error if offset is not valid under HGVS rules."""
    crossmap = NonCoding(_exons, length=75)
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=0, region='u'))
        assert e.value.args[0] == "Offset 0 at upstream boundary should be negative."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=1, region='u'))
        assert e.value.args[0] == "Offset 1 at upstream boundary should be negative."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=-6, region='u'))
        assert e.value.args[0] == "Offset -6 exceeds upstream boundary."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=-1, region=''))
        assert e.value.args[0] == "Offset -1 at the first exon should be in the upstream region."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=1, region=''))
        assert e.value.args[0] == "Offset 1 should be at a locus end."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=10, offset=1, region=''))
        assert e.value.args[0] == "Offset 1 should be at a locus end."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=10, offset=-11, region=''))
        assert e.value.args[0] == "Offset -11 exceeds intron length."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=22, offset=-1, region=''))
        assert e.value.args[0] == "Offset -1 should be at a locus start."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=22, offset=1, region=''))
        assert e.value.args[0] == "Offset 1 at the first exon should be in the downstream region."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=22, offset=0, region='d'))
        assert e.value.args[0] == "Offset 0 at downstream boundary should be positive."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=22, offset=-1, region='d'))
        assert e.value.args[0] == "Offset -1 at downstream boundary should be positive."


def test_NonCoding_invalid_offset_inverted():
    """Raise error if offset is not valid under HGVS rules."""
    crossmap = NonCoding(_exons, length=75, inverted=True)
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=0, region='u'))
        assert e.value.args[0] == "Offset 0 at upstream boundary should be negative."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=1, region='u'))
        assert e.value.args[0] == "Offset 1 at upstream boundary should be negative."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=-5, region='u'))
        assert e.value.args[0] == "Offset -5 exceeds upstream boundary."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=-1, region=''))
        assert e.value.args[0] == "Offset -1 at the first exon should be in the upstream region."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=1, offset=1, region=''))
        assert e.value.args[0] == "Offset 1 should be at a locus end."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=13, offset=-1, region=''))
        assert e.value.args[0] == "Offset -1 should be at a locus end."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=13, offset=11, region=''))
        assert e.value.args[0] == "Offset 11 exceeds intron length."
    with pytest.raises(IndexError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=22, offset=-1, region=''))
        assert e.value.args[0] == "Offset -1 should be at a locus start."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=22, offset=1, region=''))
        assert e.value.args[0] == "Offset 1 at the last exon on the reverse complement should be in the downstream region."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=22, offset=0, region='d'))
        assert e.value.args[0] == "Offset 0 at downstream boundary should be positive."
    with pytest.raises(ValueError) as e:
        crossmap.noncoding_to_coordinate(NonCodingPoint(position=22, offset=-1, region='d'))
        assert e.value.args[0] == "Offset -1 at downstream boundary should be positive."


def test_CodingPoint_invalid_initialization():
    """Raise error with invalid initialization."""
    with pytest.raises(ValueError) as e:
        CodingPoint(position=0, offset=0, region='-')
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        CodingPoint(position=0, offset=0, region='*')
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        CodingPoint(position=0, offset=0, region='')
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        CodingPoint(position=-1, offset=0, region='')
    assert str(e.value) == "Position -1 must be a positive integer."
    with pytest.raises(ValueError) as e:
        CodingPoint(position=1, offset=None, region='')
    assert str(e.value) == "Value must be an integer."
    with pytest.raises(ValueError) as e:
        CodingPoint(position=2, offset=1, region='upstream')
    assert str(e.value) == "Region must be a string in ['', 'u', 'd', '-', '*']."


def test_Coding_invalid():
    """Raise ValueError if coding is invalid."""

    with pytest.raises(ValueError) as e:
        Coding([(20, 20)], (20, 20))
    assert str(e.value) == "Locus start 20 must be smaller than locus end 20."
    with pytest.raises(ValueError) as e:
        Coding([(10, 20)], (9,15))
    assert str(e.value) == "Coordinate 9 of CDS (9, 15) is not within any exon."
    with pytest.raises(ValueError) as e:
        Coding([(10, 20)], (10,21))
    assert str(e.value) == "Coordinate 21 of CDS (10, 21) is not within any exon."
    with pytest.raises(ValueError) as e:
        Coding([(10, 20)], (15, 10))
    assert str(e.value) == "Locus start 15 must be smaller than locus end 10."
    with pytest.raises(ValueError) as e:
        Coding([], None)
    assert str(e.value) == "Locations must be a non-empty list of tuples."

    # Reverse orientation
    with pytest.raises(ValueError) as e:
        Coding([(20, 20)], (20, 20), inverted=True)
    assert str(e.value) == "Locus start 20 must be smaller than locus end 20."
    with pytest.raises(ValueError) as e:
        Coding([(10, 20)], (9,15), inverted=True)
    assert str(e.value) == "Coordinate 9 of CDS (9, 15) is not within any exon."
    with pytest.raises(ValueError) as e:
        Coding([(10, 20)], (10,21), inverted=True)
    assert str(e.value) == "Coordinate 21 of CDS (10, 21) is not within any exon."
    with pytest.raises(ValueError) as e:
        Coding([(10, 20)], (15, 10), inverted=True)
    assert str(e.value) == "Locus start 15 must be smaller than locus end 10."
    with pytest.raises(ValueError) as e:
        Coding([], None, inverted=True)
    assert str(e.value) == "Locations must be a non-empty list of tuples."


def test_Coding_invalid_with_length():
    """Raise ValueError if coordinate is out of bounds."""
    with pytest.raises(ValueError) as e:
        Coding(_exons, _cds, length=70)
    assert str(e.value) == "Value 72 must be within the bounds of the reference length 70."
    # Reverse orientation
    with pytest.raises(ValueError) as e:
        Coding(_exons, _cds, length=70, inverted=True)
    assert str(e.value) == "Value 72 must be within the bounds of the reference length 70."


def test_Coding():
    """Forward oriented coding transcript."""
    crossmap = Coding(_exons, _cds)

    # Boundary between upstream and 5' UTR.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(4),
        crossmap.coding_to_coordinate,
        CodingPoint(position=11, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(5),
        crossmap.coding_to_coordinate,
        CodingPoint(position=11, offset=0, region='-'),
    )

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

    # Boundary between 3' and downstream.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(71),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=0, region='*'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(72),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=1, region='d'),
    )


def test_Coding_with_length():
    """Raise ValueError if coordinate is out of bounds."""
    crossmap = Coding(_exons, _cds, length=75)
    # Boundary between upstream and 5' UTR.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(4),
        crossmap.coding_to_coordinate,
        CodingPoint(position=11, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(5),
        crossmap.coding_to_coordinate,
        CodingPoint(position=11, offset=0, region='-'),
    )

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

    # Boundary between 3' and downstream.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(71),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=0, region='*'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(72),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=1, region='d'),
    )

    # Boundary between downstream and sequence end.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(74),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=3, region='d'),
    )
    with pytest.raises(ValueError) as e:
        crossmap.coordinate_to_coding(Coord(75))
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=5, offset=4, region='d'))


def test_Coding_inverted():
    """Reverse oriented coding transcript."""
    crossmap = Coding(_exons, _cds, inverted=True)

    # Boundary between upstream and 5' UTR.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(72),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(71),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=0, region='-'),
    )

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

    # Boundary between 3' and downstream.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(5),
        crossmap.coding_to_coordinate,
        CodingPoint(position=11, offset=0, region='*'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(4),
        crossmap.coding_to_coordinate,
        CodingPoint(position=11, offset=1, region='d'),
    )


def test_Coding_inverted_with_length():
    """Reverse oriented coding transcript."""
    crossmap = Coding(_exons, _cds, length=75, inverted=True)

    # Boundary between upstream and sequence end.
    with pytest.raises(ValueError) as e:
        crossmap.coordinate_to_coding(Coord(75))
    assert str(e.value) == "Value 75 must be within the bounds of the reference length 75."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=5, offset=-4, region='u'))
    assert str(e.value) == "Offset -4 exceeds upstream region."
    invariant(
        crossmap.coordinate_to_coding,
        Coord(74),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=-3, region='u'),
    )

    # Boundary between upstream and 5' UTR.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(72),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=-1, region='u'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(71),
        crossmap.coding_to_coordinate,
        CodingPoint(position=5, offset=0, region='-'),
    )

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

    # Boundary between 3' and downstream.
    invariant(
        crossmap.coordinate_to_coding,
        Coord(5),
        crossmap.coding_to_coordinate,
        CodingPoint(position=11, offset=0, region='*'),
    )
    invariant(
        crossmap.coordinate_to_coding,
        Coord(4),
        crossmap.coding_to_coordinate,
        CodingPoint(position=11, offset=1, region='d'),
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
        ],
    )
    degenerate_equal(
        crossmap.coding_to_coordinate,
        Coord(9),
        [
            CodingPoint(position=1, offset=1, region='d'),
            CodingPoint(position=1, offset=0, region='*'),
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


def test_Coding_invalid_position():
    """Raise error if position in coding point is invalid."""
    crossmap = Coding(_exons, _cds)

    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=0, offset=1, region='u'))
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=12, offset=-1, region='u'))
    assert str(e.value) == "Position 12 is not in upstream boundary."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=13, offset=1, region='-'))
    assert str(e.value) == "Position 13 exceeds - region."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=-1, offset=0, region='-'))
    assert str(e.value) == "Position -1 must be a positive integer."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=0, offset=0, region=''))
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=7, offset=0, region=''))
    assert str(e.value) == "Position 7 exceeds coding region."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=6, offset=1, region='*'))
    assert str(e.value) == "Position 6 exceeds * region."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=None, offset=0, region='*'))
    assert str(e.value) == "Value must be an integer."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=6, offset=0, region='d'))
    assert str(e.value) == "Position 6 is not in downstream boundary."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=1000, offset=2, region='d'))
    assert str(e.value) == "Position 1000 is not in downstream boundary."


def test_Coding_inverted_invalid_position_inverted():
    """Raise error if position in coding point is invalid for inverted coding."""
    crossmap = Coding(_exons, _cds, inverted=True)

    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=0, offset=1, region='u'))
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=12, offset=-1, region='u'))
    assert str(e.value) == "Position 12 is not in upstream boundary."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=13, offset=1, region='-'))
    assert str(e.value) == "Position 13 exceeds - region."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=-1, offset=0, region='-'))
    assert str(e.value) == "Position -1 must be a positive integer."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=0, offset=0, region=''))
    assert str(e.value) == "Position 0 must be a positive integer."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=7, offset=0, region=''))
    assert str(e.value) == "Position 7 exceeds coding region."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=13, offset=1, region='*'))
    assert str(e.value) == "Position 13 exceeds * region."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=None, offset=0, region='*'))
    assert str(e.value) == "Value must be an integer."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=6, offset=0, region='d'))
    assert str(e.value) == "Position 6 is not in downstream boundary."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=1000, offset=2, region='d'))
    assert str(e.value) == "Position 1000 is not in downstream boundary."


def test_Coding_invalid_offset():
    """Raise error if offset in coding point is invalid."""
    crossmap = Coding(_exons, _cds, length=75)

    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=1, offset=-6, region='u'))
    assert str(e.value) == "Offset -6 exceeds upstream region."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=1, offset=1, region='-'))
    assert str(e.value) == "Position 1 is not at a locus boundary."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=1, offset=-1, region='-'))
    assert str(e.value) == "Position 1 is not at a locus boundary."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=1, offset=-1, region=''))
    assert str(e.value) == "Position 1 is not at a locus boundary."
    with pytest.raises(IndexError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=3, offset=6, region=''))
    assert str(e.value) == "Offset 6 exceeds intron length."
    with pytest.raises(IndexError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=4, offset=12, region='*'))
    assert str(e.value) == "Offset 12 should be at a locus end."
    with pytest.raises(IndexError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=4, offset=-50, region='*'))
    assert str(e.value) == "Offset -50 exceeds intron length."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=5, offset=10, region='d'))
    assert str(e.value) == "Offset 10 exceeds downstream region."


def test_Coding_invalid_offset_inverted():
    """Raise error if offset in coding point is invalid."""
    crossmap = Coding(_exons, _cds, length=75, inverted=True)

    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=1, offset=-6, region='u'))
    assert str(e.value) == "Offset -6 exceeds upstream region."
    with pytest.raises(IndexError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=1, offset=1, region='-'))
    assert str(e.value) == "Offset 1 should be at a locus end."
    with pytest.raises(IndexError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=2, offset=-1, region='-'))
    assert str(e.value) == "Offset -1 should be at a locus start."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=1, offset=-1, region=''))
    assert str(e.value) == "Position 1 is not at a locus boundary."
    with pytest.raises(IndexError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=3, offset=6, region=''))
    assert str(e.value) == "Offset 6 exceeds intron length."
    with pytest.raises(IndexError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=4, offset=12, region='*'))
    assert str(e.value) == "Offset 12 exceeds intron length."
    with pytest.raises(IndexError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=4, offset=-50, region='*'))
    assert str(e.value) == "Offset -50 exceeds intron length."
    with pytest.raises(ValueError) as e:
        crossmap.coding_to_coordinate(CodingPoint(position=11, offset=10, region='d'))
    assert str(e.value) == "Offset 10 exceeds downstream region."


def test_Coding_protein():
    """Protein positions."""
    crossmap = Coding(_exons, _cds)

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
