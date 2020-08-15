from mutalyzer_crossmapper import Coding, Genomic, NonCoding

from helper import degenerate_equal, invariant

_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_cds = (32, 43)


def test_Genomic():
    """Genomic positions are coordinates incremented by one."""
    crossmap = Genomic()

    invariant(
        crossmap.coordinate_to_genomic, 0, crossmap.genomic_to_coordinate, 1)
    invariant(
        crossmap.coordinate_to_genomic, 98, crossmap.genomic_to_coordinate, 99)


def test_NonCoding():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons)

    # Boundary between upstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding, 4,
        crossmap.noncoding_to_coordinate, (1, -1, -1))
    invariant(
        crossmap.coordinate_to_noncoding, 5,
        crossmap.noncoding_to_coordinate, (1, 0, 0))

    # Boundary between downstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding, 71,
        crossmap.noncoding_to_coordinate, (22, 0, 0))
    invariant(
        crossmap.coordinate_to_noncoding, 72,
        crossmap.noncoding_to_coordinate, (22, 1, 1))


def test_NonCoding_inverted():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons, inverted=True)

    # Boundary between upstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding, 72,
        crossmap.noncoding_to_coordinate, (1, -1, -1))
    invariant(
        crossmap.coordinate_to_noncoding, 71,
        crossmap.noncoding_to_coordinate, (1, 0, 0))

    # Boundary between downstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding, 5,
        crossmap.noncoding_to_coordinate, (22, 0, 0))
    invariant(
        crossmap.coordinate_to_noncoding, 4,
        crossmap.noncoding_to_coordinate, (22, 1, 1))


def test_NonCoding_degenerate():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons)

    # Boundary between upstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate, 4,
        [(1, -1, -1), (-1, 0, -1)])

    # Boundary between downstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate, 72,
        [(22, 1, 1), (23, 0, 1)])


def test_NonCoding_inverted_degenerate():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons, inverted=True)

    # Boundary between upstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate, 72,
        [(1, -1, -1), (-1, 0, -1)])

    # Boundary between downstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate, 4,
        [(22, 1, 1), (23, 0, 1)])


def test_Coding():
    """Forward oriented coding transcript."""
    crossmap = Coding(_exons, _cds)

    # Boundary between 5' and CDS.
    invariant(
        crossmap.coordinate_to_coding, 31,
        crossmap.coding_to_coordinate, (-1, 0, -1, 0))
    invariant(
        crossmap.coordinate_to_coding, 32,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))

    # Boundary between CDS and 3'.
    invariant(
        crossmap.coordinate_to_coding, 42,
        crossmap.coding_to_coordinate, (6, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (1, 0, 1, 0))


def test_Coding_inverted():
    """Reverse oriented coding transcript."""
    crossmap = Coding(_exons, _cds, True)

    # Boundary between 5' and CDS.
    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (-1, 0, -1, 0))
    invariant(
        crossmap.coordinate_to_coding, 42,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))

    # Boundary between CDS and 3'.
    invariant(
        crossmap.coordinate_to_coding, 32,
        crossmap.coding_to_coordinate, (6, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 31,
        crossmap.coding_to_coordinate, (1, 0, 1, 0))


def test_Coding_regions():
    """The CDS can start or end on a region boundary."""
    crossmap = Coding([(10, 21), (30, 40), (49, 60)], (30, 40))

    # Upstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding, 25,
        crossmap.coding_to_coordinate, (-1, 5, -1, 0))
    invariant(
        crossmap.coordinate_to_coding, 26,
        crossmap.coding_to_coordinate, (1, -4, 0, 0))

    # Downstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding, 44,
        crossmap.coding_to_coordinate, (10, 5, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 45,
        crossmap.coding_to_coordinate, (1, -4, 1, 0))


def test_Coding_regions_inverted():
    """The CDS can start or end on a region boundary."""
    crossmap = Coding([(10, 21), (30, 40), (49, 60)], (30, 40), True)

    # Upstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding, 44,
        crossmap.coding_to_coordinate, (-1, 5, -1, 0))
    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (1, -4, 0, 0))

    # Downstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding, 25,
        crossmap.coding_to_coordinate, (10, 5, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 24,
        crossmap.coding_to_coordinate, (1, -4, 1, 0))


def test_Coding_no_utr5():
    """A 5' UTR may be missing."""
    crossmap = Coding([(10, 20)], (10, 15))

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (1, -1, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))


def test_Coding_no_utr5_inverted():
    """A 5' UTR may be missing."""
    crossmap = Coding([(10, 20)], (15, 20), True)

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (1, -1, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))


def test_Coding_no_utr3():
    """A 3' UTR may be missing."""
    crossmap = Coding([(10, 20)], (15, 20))

    # Direct transition from CDS to downstream.
    invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (5, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (5, 1, 0, 1))


def test_Coding_no_utr3_inverted():
    """A 3' UTR may be missing."""
    crossmap = Coding([(10, 20)], (10, 15), True)

    # Direct transition from CDS to downstream.
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (5, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (5, 1, 0, 1))


def test_Coding_small_utr5():
    """A 5' UTR may be of lenght one."""
    crossmap = Coding([(10, 20)], (11, 15))

    # Transition from upstream to 5' UTR to CDS.
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (-1, -1, -1, -1))
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (-1, 0, -1, 0))
    invariant(
        crossmap.coordinate_to_coding, 11,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))


def test_Coding_small_utr5_inverted():
    """A 5' UTR may be of lenght one."""
    crossmap = Coding([(10, 20)], (15, 19), True)

    # Transition from upstream to 5' UTR to CDS.
    invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (-1, -1, -1, -1))
    invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (-1, 0, -1, 0))
    invariant(
        crossmap.coordinate_to_coding, 18,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))


def test_Coding_small_utr3():
    """A 5' UTR may be of lenght one."""
    crossmap = Coding([(10, 20)], (15, 19))

    # Transition from CDS to 3' UTR to downstream.
    invariant(
        crossmap.coordinate_to_coding, 18,
        crossmap.coding_to_coordinate, (4, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (1, 0, 1, 0))
    invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (1, 1, 1, 1))


def test_Coding_small_utr3_inverted():
    """A 5' UTR may be of lenght one."""
    crossmap = Coding([(10, 20)], (11, 15), True)

    # Transition from CDS to 3' UTR to downstream.
    invariant(
        crossmap.coordinate_to_coding, 11,
        crossmap.coding_to_coordinate, (4, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (1, 0, 1, 0))
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (1, 1, 1, 1))


def test_Coding_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    crossmap = Coding([(10, 20)], (11, 19))

    degenerate_equal(
        crossmap.coding_to_coordinate, 9,
        [(-1, -1, -1, -1), (-2, 0, -1, -1), (1, -2, 0, -1), (1, -10, 1, -1)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 20,
        [(1, 1, 1, 1), (2, 0, 1, 1), (8, 2, 0, 1), (-1, 10, -1, 1)])


def test_Coding_inverted_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    crossmap = Coding([(10, 20)], (11, 19), True)

    degenerate_equal(
        crossmap.coding_to_coordinate, 20,
        [(-1, -1, -1, -1), (-2, 0, -1, -1), (1, -2, 0, -1), (1, -10, 1, -1)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 9,
        [(1, 1, 1, 1), (2, 0, 1, 1), (8, 2, 0, 1), (-1, 10, -1, 1)])


def test_Coding_degenerate_return():
    """Degenerate upstream and downstream positions may be returned."""
    crossmap = Coding([(10, 20)], (11, 19))

    assert crossmap.coordinate_to_coding(9, True) == (-2, 0, -1, -1)
    assert crossmap.coordinate_to_coding(20, True) == (2, 0, 1, 1)


def test_Coding_inverted_degenerate_return():
    """Degenerate upstream and downstream positions may be returned."""
    crossmap = Coding([(10, 20)], (11, 19), True)

    assert crossmap.coordinate_to_coding(20, True) == (-2, 0, -1, -1)
    assert crossmap.coordinate_to_coding(9, True) == (2, 0, 1, 1)


def test_Coding_degenerate_no_return():
    """Degenerate internal positions do not exist."""
    crossmap = Coding([(10, 20), (30, 40)], (10, 40))

    assert (crossmap.coordinate_to_coding(25) ==
            crossmap.coordinate_to_coding(25, True))


def test_Coding_inverted_degenerate_no_return():
    """Degenerate internal positions do not exist."""
    crossmap = Coding([(10, 20), (30, 40)], (10, 40), True)

    assert (crossmap.coordinate_to_coding(25) ==
            crossmap.coordinate_to_coding(25, True))


def test_Coding_no_utr_degenerate():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11))

    degenerate_equal(
        crossmap.coding_to_coordinate, 9,
        [(1, -1, 0, -1), (-1, 0, -1, -1), (1, -2, 1, -1)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 11,
        [(1, 1, 0, 1), (1, 0, 1, 1), (-1, 2, -1, 1)])


def test_Coding_inverted_no_utr_degenerate():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11), True)

    degenerate_equal(
        crossmap.coding_to_coordinate, 11,
        [(1, -1, 0, -1), (-1, 0, -1, -1), (1, -2, 1, -1)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 9,
        [(1, 1, 0, 1), (1, 0, 1, 1), (-1, 2, -1, 1)])


def test_Coding_no_utr_degenerate_return():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11))

    assert crossmap.coordinate_to_coding(8, True) == (-2, 0, -1, -2)
    assert crossmap.coordinate_to_coding(9, True) == (-1, 0, -1, -1)
    assert crossmap.coordinate_to_coding(11, True) == (1, 0, 1, 1)
    assert crossmap.coordinate_to_coding(12, True) == (2, 0, 1, 2)


def test_Coding_inverted_no_utr_degenerate_return():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11), True)

    assert crossmap.coordinate_to_coding(11, True) == (-1, 0, -1, -1)
    assert crossmap.coordinate_to_coding(9, True) == (1, 0, 1, 1)


def test_Coding_protein():
    """Protein positions."""
    crossmap = Coding(_exons, _cds)

    # Boundary between 5' UTR and CDS.
    invariant(
        crossmap.coordinate_to_protein, 31,
        crossmap.protein_to_coordinate, (-1, 3, 0, -1, 0))
    invariant(
        crossmap.coordinate_to_protein, 32,
        crossmap.protein_to_coordinate, (1, 1, 0, 0, 0))

    # Intron boundary.
    invariant(
        crossmap.coordinate_to_protein, 34,
        crossmap.protein_to_coordinate, (1, 3, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_protein, 35,
        crossmap.protein_to_coordinate, (1, 3, 1, 0, 0))

    # Boundary between CDS and 3' UTR.
    invariant(
        crossmap.coordinate_to_protein, 42,
        crossmap.protein_to_coordinate, (2, 3, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_protein, 43,
        crossmap.protein_to_coordinate, (1, 1, 0, 1, 0))
