from mutalyzer_crossmapper import Crossmap

from helper import degenerate_equal, invariant

_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_cds = (32, 43)


def test_Crossmap_genomic():
    """Genomic positions are coordinates incremented by one."""
    crossmap = Crossmap()

    invariant(
        crossmap.coordinate_to_genomic, 0, crossmap.genomic_to_coordinate, 1)
    invariant(
        crossmap.coordinate_to_genomic, 98, crossmap.genomic_to_coordinate, 99)


#def test_Locus_degenerate_return():
#    """Degenerate position can be retured."""
#    locus = Locus((10, 20))
#
#    assert locus.to_position(9, True) == (-1, 0)
#    assert locus.to_position(20, True) == (10, 0)


#def test_Locus_inverted_degenerate_return():
#    """Degenerate position can be retured."""
#    locus = Locus((10, 20), True)
#
#    assert locus.to_position(20, True) == (-1, 0)
#    assert locus.to_position(9, True) == (11, 0)


#def test_MultiLocus_degenerate_return():
#    """Degenerate upstream and downstream positions may be returned."""
#    multi_locus = MultiLocus(_locations)
#
#    assert multi_locus.to_position(4, True) == (-1, 0, -1)
#    assert multi_locus.to_position(72, True) == (23, 0, 1)
#
#
#def test_MultiLocus_inverted_degenerate_return():
#    """Degenerate upstream and downstream positions may be returned."""
#    multi_locus = MultiLocus(_locations, True)
#
#    assert multi_locus.to_position(72, True) == (-1, 0, -1)
#    assert multi_locus.to_position(4, True) == (23, 0, 1)
#
#
#def test_MultiLocus_degenerate_no_return():
#    """Degenerate internal positions do not exist."""
#    multi_locus = MultiLocus(_locations)
#
#    assert multi_locus.to_position(10, True) == multi_locus.to_position(10)
#
#
#def test_MultiLocus_inverted_degenerate_no_return():
#    """Degenerate internal positions do not exist."""
#    multi_locus = MultiLocus(_locations, True)
#
#    assert multi_locus.to_position(10, True) == multi_locus.to_position(10)


def test_Crossmap_noncoding():
    """Forward oriented noncoding transcript."""
    crossmap = Crossmap(_exons)

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


def test_Crossmap_noncoding_inverted():
    """Forward oriented noncoding transcript."""
    crossmap = Crossmap(_exons, inverted=True)

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


def test_Crossmap_noncoding_degenerate():
    """Forward oriented noncoding transcript."""
    crossmap = Crossmap(_exons)

    # Boundary between upstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate, 4,
        [(1, -1, -1), (-1, 0, -1)])

    # Boundary between downstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate, 72,
        [(22, 1, 1), (23, 0, 1)])


def test_Crossmap_noncoding_inverted_degenerate():
    """Forward oriented noncoding transcript."""
    crossmap = Crossmap(_exons, inverted=True)

    # Boundary between upstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate, 72,
        [(1, -1, -1), (-1, 0, -1)])

    # Boundary between downstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate, 4,
        [(22, 1, 1), (23, 0, 1)])


def test_Crossmap_coding():
    """Forward oriented coding transcript."""
    crossmap = Crossmap(_exons, _cds)

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


def test_Crossmap_coding_inverted():
    """Reverse oriented coding transcript."""
    crossmap = Crossmap(_exons, _cds, True)

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


def test_Crossmap_coding_regions():
    """The CDS can start or end on a region boundary."""
    crossmap = Crossmap([(10, 21), (30, 40), (49, 60)], (30, 40))

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


def test_Crossmap_coding_regions_inverted():
    """The CDS can start or end on a region boundary."""
    crossmap = Crossmap([(10, 21), (30, 40), (49, 60)], (30, 40), True)

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


def test_Crossmap_coding_no_utr5():
    """A 5' UTR may be missing."""
    crossmap = Crossmap([(10, 20)], (10, 15))

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (1, -1, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))


def test_Crossmap_coding_no_utr5_inverted():
    """A 5' UTR may be missing."""
    crossmap = Crossmap([(10, 20)], (15, 20), True)

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (1, -1, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))


def test_Crossmap_coding_no_utr3():
    """A 3' UTR may be missing."""
    crossmap = Crossmap([(10, 20)], (15, 20))

    # Direct transition from CDS to downstream.
    invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (5, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (5, 1, 0, 1))


def test_Crossmap_coding_no_utr3_inverted():
    """A 3' UTR may be missing."""
    crossmap = Crossmap([(10, 20)], (10, 15), True)

    # Direct transition from CDS to downstream.
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (5, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (5, 1, 0, 1))


def test_Crossmap_coding_small_utr5():
    """A 5' UTR may be of lenght one."""
    crossmap = Crossmap([(10, 20)], (11, 15))

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


def test_Crossmap_coding_small_utr5_inverted():
    """A 5' UTR may be of lenght one."""
    crossmap = Crossmap([(10, 20)], (15, 19), True)

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


def test_Crossmap_coding_small_utr3():
    """A 5' UTR may be of lenght one."""
    crossmap = Crossmap([(10, 20)], (15, 19))

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


def test_Crossmap_coding_small_utr3_inverted():
    """A 5' UTR may be of lenght one."""
    crossmap = Crossmap([(10, 20)], (11, 15), True)

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


def test_Crossmap_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    crossmap = Crossmap([(10, 20)], (11, 19))

    degenerate_equal(
        crossmap.coding_to_coordinate, 9,
        [(-1, -1, -1, -1), (-2, 0, -1, -1), (1, -2, 0, -1), (1, -10, 1, -1)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 20,
        [(1, 1, 1, 1), (2, 0, 1, 1), (8, 2, 0, 1), (-1, 10, -1, 1)])


def test_Crossmap_inverted_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    crossmap = Crossmap([(10, 20)], (11, 19), True)

    degenerate_equal(
        crossmap.coding_to_coordinate, 20,
        [(-1, -1, -1, -1), (-2, 0, -1, -1), (1, -2, 0, -1), (1, -10, 1, -1)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 9,
        [(1, 1, 1, 1), (2, 0, 1, 1), (8, 2, 0, 1), (-1, 10, -1, 1)])


#def test_Crossmap_degenerate_return():
#    """Degenerate upstream and downstream positions may be returned."""
#    crossmap = Crossmap([(10, 20)], (11, 19))
#
#    assert crossmap.coordinate_to_coding(9, True) == (-2, 0, -1, -1)
#    assert crossmap.coordinate_to_coding(20, True) == (2, 0, 1, 1)


#def test_Crossmap_inverted_degenerate_return():
#    """Degenerate upstream and downstream positions may be returned."""
#    crossmap = Crossmap([(10, 20)], (11, 19), True)
#
#    assert crossmap.coordinate_to_coding(20, True) == (-2, 0, -1, -1)
#    assert crossmap.coordinate_to_coding(9, True) == (2, 0, 1, 1)


#def test_Crossmap_degenerate_no_return():
#    """Degenerate internal positions do not exist."""
#    crossmap = Crossmap([(10, 20), (30, 40)], (10, 40))
#
#    assert (crossmap.coordinate_to_coding(25) ==
#            crossmap.coordinate_to_coding(25, True))


#def test_Crossmap_inverted_degenerate_no_return():
#    """Degenerate internal positions do not exist."""
#    crossmap = Crossmap([(10, 20), (30, 40)], (10, 40), True)
#
#    assert (crossmap.coordinate_to_coding(25) ==
#            crossmap.coordinate_to_coding(25, True))


def test_Crossmap_no_utr_degenerate():
    """UTRs may be missing."""
    crossmap = Crossmap([(10, 11)], (10, 11))

    degenerate_equal(
        crossmap.coding_to_coordinate, 9,
        [(1, -1, 0, -1), (-1, 0, -1, -1), (1, -2, 1, -1)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 11,
        [(1, 1, 0, 1), (1, 0, 1, 1), (-1, 2, -1, 1)])


def test_Crossmap_inverted_no_utr_degenerate():
    """UTRs may be missing."""
    crossmap = Crossmap([(10, 11)], (10, 11), True)

    degenerate_equal(
        crossmap.coding_to_coordinate, 11,
        [(1, -1, 0, -1), (-1, 0, -1, -1), (1, -2, 1, -1)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 9,
        [(1, 1, 0, 1), (1, 0, 1, 1), (-1, 2, -1, 1)])


#def test_Crossmap_no_utr_degenerate_return():
#    """UTRs may be missing."""
#    crossmap = Crossmap([(10, 11)], (10, 11))
#
#    assert crossmap.coordinate_to_coding(9, True) == (-1, 0, -1, 0)
#    assert crossmap.coordinate_to_coding(11, True) == (2, 0, 1, 0)


#def test_Crossmap_inverted_no_utr_degenerate_return():
#    """UTRs may be missing."""
#    crossmap = Crossmap([(10, 11)], (10, 11), True)
#
#    assert crossmap.coordinate_to_coding(11, True) == (-1, 0, -1, 0)
#    assert crossmap.coordinate_to_coding(9, True) == (2, 0, 1, 0)


#def test_Crossmap_protein():
#    crossmap = Crossmap(_exons, _cds)
#
#    invariant(
#        crossmap.coordinate_to_protein, 31,
#        crossmap.protein_to_coordinate, (-1, 3, 0, 0))
#    invariant(
#        crossmap.coordinate_to_protein, 32,
#        crossmap.protein_to_coordinate, (1, 1, 0, 1))
#    invariant(
#        crossmap.coordinate_to_protein, 34,
#        crossmap.protein_to_coordinate, (1, 3, 0, 1))
#    invariant(
#        crossmap.coordinate_to_protein, 36,
#        crossmap.protein_to_coordinate, (1, 3, 2, 1))
#    invariant(
#        crossmap.coordinate_to_protein, 40,
#        crossmap.protein_to_coordinate, (2, 1, 0, 1))
#    invariant(
#        crossmap.coordinate_to_protein, 42,
#        crossmap.protein_to_coordinate, (2, 3, 0, 1))
#    invariant(
#        crossmap.coordinate_to_protein, 43,
#        crossmap.protein_to_coordinate, (1, 1, 0, 2))
