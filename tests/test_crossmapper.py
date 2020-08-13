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


def test_Crossmap_coding():
    """Forward oriented Crossmap."""
    crossmap = Crossmap(_exons, _cds)

    # Boundary between 5' and CDS.
    invariant(
        crossmap.coordinate_to_coding, 31,
        crossmap.coding_to_coordinate, (-1, 0, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 32,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))

    # Boundary between CDS and 3'.
    invariant(
        crossmap.coordinate_to_coding, 42,
        crossmap.coding_to_coordinate, (6, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (1, 0, 0, 1))


def test_Crossmap_coding_inverted():
    """Reverse oriented Crossmap."""
    crossmap = Crossmap(_exons, _cds, True)

    # Boundary between 5' and CDS.
    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (-1, 0, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 42,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))

    # Boundary between CDS and 3'.
    invariant(
        crossmap.coordinate_to_coding, 32,
        crossmap.coding_to_coordinate, (6, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 31,
        crossmap.coding_to_coordinate, (1, 0, 0, 1))


def test_Crossmap_coding_regions():
    """The CDS can start or end on a region boundary."""
    crossmap = Crossmap([(10, 21), (30, 40), (49, 60)], (30, 40))

    # Upstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding, 25,
        crossmap.coding_to_coordinate, (-1, 5, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 26,
        crossmap.coding_to_coordinate, (1, -4, 0, 0))

    # Downstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding, 44,
        crossmap.coding_to_coordinate, (10, 5, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 45,
        crossmap.coding_to_coordinate, (1, -4, 0, 1))


def test_Crossmap_coding_regions_inverted():
    """The CDS can start or end on a region boundary."""
    crossmap = Crossmap([(10, 21), (30, 40), (49, 60)], (30, 40), True)

    # Upstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding, 44,
        crossmap.coding_to_coordinate, (-1, 5, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (1, -4, 0, 0))

    # Downstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding, 25,
        crossmap.coding_to_coordinate, (10, 5, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 24,
        crossmap.coding_to_coordinate, (1, -4, 0, 1))


def test_Crossmap_coding_no_utr5():
    """A 5' UTR may be missing."""
    crossmap = Crossmap([(10, 20)], (10, 15))

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (1, -1, -1, 0))
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))


def test_Crossmap_coding_no_utr5_inverted():
    """A 5' UTR may be missing."""
    crossmap = Crossmap([(10, 20)], (15, 20), True)

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (1, -1, -1, 0))
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
        crossmap.coding_to_coordinate, (5, 1, 1, 0))


def test_Crossmap_coding_no_utr3_inverted():
    """A 3' UTR may be missing."""
    crossmap = Crossmap([(10, 20)], (10, 15), True)

    # Direct transition from CDS to downstream.
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (5, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (5, 1, 1, 0))


def test_Crossmap_coding_small_utr5():
    """A 5' UTR may be of lenght one."""
    crossmap = Crossmap([(10, 20)], (11, 15))

    # Transition from upstream to 5' UTR to CDS.
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (-1, -1, -1, -1))
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (-1, 0, 0, -1))
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
        crossmap.coding_to_coordinate, (-1, 0, 0, -1))
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
        crossmap.coding_to_coordinate, (1, 0, 0, 1))
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
        crossmap.coding_to_coordinate, (1, 0, 0, 1))
    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (1, 1, 1, 1))


#def test_Crossmap_degenerate():
#    """Degenerate upstream and downstream positions are silently corrected."""
#    crossmap = Crossmap([(10, 20)], (11, 19))
#
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 9,
#        [(-1, -1, -1, -1), (-2, 0, -1, -1), (-2, 0, -2, 0), (1, -2, -1, 0)])
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 20,
#        [(1, 1, 1, 1), (2, 0, 1, 1), (10, 0, 2, 0), (8, 2, 2, 0)])
#
#
#def test_Crossmap_inverted_degenerate():
#    """Degenerate upstream and downstream positions are silently corrected."""
#    crossmap = Crossmap([(10, 20)], (11, 19), True)
#
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 20,
#        [(-1, -1, -1, -1), (-2, 0, -1, -1), (-2, 0, -2, 0), (1, -2, -1, 0)])
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 9,
#        [(1, 1, 1, 1), (2, 0, 1, 1), (10, 0, 2, 0), (8, 2, 2, 0)])
#
#
#def test_Crossmap_degenerate_return():
#    """Degenerate upstream and downstream positions may be returned."""
#    crossmap = Crossmap([(10, 20)], (11, 19))
#
#    assert crossmap.coordinate_to_coding(9, True) == (-2, 0, -1, -1)
#    assert crossmap.coordinate_to_coding(20, True) == (2, 0, 1, 1)
#
#
#def test_Crossmap_inverted_degenerate_return():
#    """Degenerate upstream and downstream positions may be returned."""
#    crossmap = Crossmap([(10, 20)], (11, 19), True)
#
#    assert crossmap.coordinate_to_coding(20, True) == (-2, 0, -1, -1)
#    assert crossmap.coordinate_to_coding(9, True) == (2, 0, 1, 1)
#
#
#def test_Crossmap_degenerate_no_return():
#    """Degenerate internal positions do not exist."""
#    crossmap = Crossmap([(10, 20), (30, 40)], (10, 40))
#
#    assert (crossmap.coordinate_to_coding(25) ==
#            crossmap.coordinate_to_coding(25, True))
#
#
#def test_Crossmap_inverted_degenerate_no_return():
#    """Degenerate internal positions do not exist."""
#    crossmap = Crossmap([(10, 20), (30, 40)], (10, 40), True)
#
#    assert (crossmap.coordinate_to_coding(25) ==
#            crossmap.coordinate_to_coding(25, True))
#
#
#def test_Crossmap_regions_degenerate_no_return():
#    """Degenerate positions do not exist between regions."""
#    crossmap = Crossmap([(10, 21), (30, 40), (49, 60)], (30, 40))
#
#    assert (crossmap.coordinate_to_coding(25) ==
#            crossmap.coordinate_to_coding(25, True))
#
#
#def test_Crossmap_inverted_regions_degenerate_no_return():
#    """Degenerate positions do not exist between regions."""
#    crossmap = Crossmap([(10, 21), (30, 40), (49, 60)], (30, 40), True)
#
#    assert (crossmap.coordinate_to_coding(25) ==
#            crossmap.coordinate_to_coding(25, True))
#
#
#def test_Crossmap_no_utr_degenerate():
#    """UTRs may be missing."""
#    crossmap = Crossmap([(10, 11)], (10, 11))
#
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 9, [(1, -1, -1, 0), (-1, 0, -1, 0)])
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 11, [(1, 1, 1, 0), (2, 0, 1, 0)])
#
#
#def test_Crossmap_no_utr_degenerate_return():
#    """UTRs may be missing."""
#    crossmap = Crossmap([(10, 11)], (10, 11))
#
#    assert crossmap.coordinate_to_coding(9, True) == (-1, 0, -1, 0)
#    assert crossmap.coordinate_to_coding(11, True) == (2, 0, 1, 0)
#
#
#def test_Crossmap_inverted_no_utr_degenerate():
#    """UTRs may be missing."""
#    crossmap = Crossmap([(10, 11)], (10, 11), True)
#
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 11, [(1, -1, -1, 0), (-1, 0, -1, 0)])
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 9, [(1, 1, 1, 0), (2, 0, 1, 0)])
#
#
#def test_Crossmap_inverted_no_utr_degenerate_return():
#    """UTRs may be missing."""
#    crossmap = Crossmap([(10, 11)], (10, 11), True)
#
#    assert crossmap.coordinate_to_coding(11, True) == (-1, 0, -1, 0)
#    assert crossmap.coordinate_to_coding(9, True) == (2, 0, 1, 0)
#
#
#def test_Crossmap_small_utr_degenerate():
#    """UTRs may be small."""
#    crossmap = Crossmap([(9, 12)], (10, 11))
#
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 8, [(-1, -1, -1, -1), (-2, 0, -1, -1)])
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 12, [(1, 1, 1, 1), (2, 0, 1, 1)])
#
#def test_Crossmap_small_utr_degenerate_return():
#    """UTRs may be small."""
#    crossmap = Crossmap([(9, 12)], (10, 11))
#
#    assert crossmap.coordinate_to_coding(8, True) == (-2, 0, -1, -1)
#    assert crossmap.coordinate_to_coding(12, True) == (2, 0, 1, 1)
#
#
#def test_Crossmap_inverted_small_utr_degenerate():
#    """UTRs may be small."""
#    crossmap = Crossmap([(9, 12)], (10, 11), True)
#
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 12, [(-1, -1, -1, -1), (-2, 0, -1, -1)])
#    degenerate_equal(
#        crossmap.coding_to_coordinate, 8, [(1, 1, 1, 1), (2, 0, 1, 1)])
#
#
#def test_Crossmap_inverted_small_utr_degenerate_return():
#    """UTRs may be small."""
#    crossmap = Crossmap([(9, 12)], (10, 11), True)
#
#    assert crossmap.coordinate_to_coding(12, True) == (-2, 0, -1, -1)
#    assert crossmap.coordinate_to_coding(8, True) == (2, 0, 1, 1)


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
