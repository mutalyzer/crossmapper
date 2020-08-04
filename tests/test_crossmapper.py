from mutalyzer_crossmapper import Crossmap

from helper import degenerate_equal, invariant


_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_cds = (32, 43)


def test_Crossmap_genomic():
    crossmap = Crossmap()

    invariant(
        crossmap.coordinate_to_genomic, 0, crossmap.genomic_to_coordinate, 1)
    invariant(
        crossmap.coordinate_to_genomic, 98, crossmap.genomic_to_coordinate, 99)


def test_Crossmap_coding():
    crossmap = Crossmap(_exons, _cds)

    invariant(
        crossmap.coordinate_to_coding, 31,
        crossmap.coding_to_coordinate, (-1, 0, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 32,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (1, 0, 0, 1))


def test_Crossmap_coding_inverted():
    crossmap = Crossmap(_exons, _cds, True)

    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (-1, 0, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 42,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 31,
        crossmap.coding_to_coordinate, (1, 0, 0, 1))


def test_Crossmap_coding_regions():
    crossmap = Crossmap([(10, 21), (30, 40), (49, 60)], (30, 40))

    invariant(
        crossmap.coordinate_to_coding, 24,
        crossmap.coding_to_coordinate, (-1, 4, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 25,
        crossmap.coding_to_coordinate, (-1, 5, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 26,
        crossmap.coding_to_coordinate, (1, -4, 0, 0))

    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (10, 4, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 44,
        crossmap.coding_to_coordinate, (10, 5, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 45,
        crossmap.coding_to_coordinate, (1, -4, 0, 1))


def test_Crossmap_coding_inverted_regions():
    crossmap = Crossmap([(10, 21), (30, 40), (49, 60)], (30, 40), True)

    invariant(
        crossmap.coordinate_to_coding, 24,
        crossmap.coding_to_coordinate, (1, -4, 0, 1))
    invariant(
        crossmap.coordinate_to_coding, 25,
        crossmap.coding_to_coordinate, (10, 5, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 26,
        crossmap.coding_to_coordinate, (10, 4, 0, 0))

    invariant(
        crossmap.coordinate_to_coding, 43,
        crossmap.coding_to_coordinate, (1, -4, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 44,
        crossmap.coding_to_coordinate, (-1, 5, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 45,
        crossmap.coding_to_coordinate, (-1, 4, 0, -1))


def test_Crossmap_degenerate():
    crossmap = Crossmap([(10, 20)], (11, 19))

    degenerate_equal(
        crossmap.coding_to_coordinate, 9, [
            (-1, -1, -1, -1), (-2, 0, -1, -1), (-2, 0, -2, 0), (1, -2, -1, 0)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 20,
        [(1, 1, 1, 1), (2, 0, 1, 1), (10, 0, 2, 0), (8, 2, 2, 0)])


def test_Crossmap_degenerate_return():
    pass


def test_Crossmap_inverted_degenerate():
    crossmap = Crossmap([(10, 20)], (11, 19), True)

    degenerate_equal(
        crossmap.coding_to_coordinate, 20, [
            (-1, -1, -1, -1), (-2, 0, -1, -1), (-2, 0, -2, 0), (1, -2, -1, 0)])
    degenerate_equal(
        crossmap.coding_to_coordinate, 9,
        [(1, 1, 1, 1), (2, 0, 1, 1), (10, 0, 2, 0), (8, 2, 2, 0)])


def test_Crossmap_inverted_degenerate_return():
    pass


#def test_Crossmap_degenerate_cds():
#    crossmap = Crossmap([(10, 11)], (10, 11))
#
#    degenerate_equal(
#        crossmap.coding_to_coordinate, [(1, 1, 0), (1, 0, 2)])
#    degenerate_equal(
#        crossmap.coding_to_coordinate, [(1, -1, 0), (-1, 0, 0)])
#
#
#def test_Crossmap_intron_degenerate():
#    crossmap = Crossmap([(10, 20), (30, 40)], (10, 40))
#
#    assert (crossmap.coordinate_to_coding(25) ==
#            crossmap.coordinate_to_coding(25, True))
#
#
#def test_Crossmap_utr_intron_degenerate():
#    crossmap = Crossmap([(0, 20), (30, 40), (50, 60)], (32, 38))
#
#    assert crossmap.coordinate_to_coding(28, True) == (-2, -2, 0, 1)
#    assert crossmap.coordinate_to_coding(42, True) == (2, 3, 2, 1)
#
#
#def test_Crossmap_utr_intron_inverted_degenerate():
#    crossmap = Crossmap([(0, 20), (30, 40), (50, 60)], (32, 38), True)
#
#    assert crossmap.coordinate_to_coding(42, True) == (-2, -3, 0, 1)
#    assert crossmap.coordinate_to_coding(28, True) == (2, 2, 2, 1)


def test_Crossmap_coding_no_utr5():
    crossmap = Crossmap([(10, 20)], (10, 15))

    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (1, -1, -1, 0))
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))


def test_Crossmap_coding_small_utr5():
    crossmap = Crossmap([(10, 20)], (11, 15))

    invariant(
        crossmap.coordinate_to_coding, 9,
        crossmap.coding_to_coordinate, (-1, -1, -1, -1))
    invariant(
        crossmap.coordinate_to_coding, 10,
        crossmap.coding_to_coordinate, (-1, 0, 0, -1))
    invariant(
        crossmap.coordinate_to_coding, 11,
        crossmap.coding_to_coordinate, (1, 0, 0, 0))


def test_Crossmap_coding_no_utr3():
    crossmap = Crossmap([(10, 20)], (15, 20))

    invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (5, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (5, 1, 1, 0))


def test_Crossmap_coding_small_utr3():
    crossmap = Crossmap([(10, 20)], (15, 19))

    invariant(
        crossmap.coordinate_to_coding, 18,
        crossmap.coding_to_coordinate, (4, 0, 0, 0))
    invariant(
        crossmap.coordinate_to_coding, 19,
        crossmap.coding_to_coordinate, (1, 0, 0, 1))
    invariant(
        crossmap.coordinate_to_coding, 20,
        crossmap.coding_to_coordinate, (1, 1, 1, 1))


#def test_Crossmap_no_utr_degenerate():
#    crossmap = Crossmap([(10, 20)], (10, 20))
#
#    degenerate_equal(
#        crossmap.coding_to_coordinate, [(1, -1, 1), (-1, 0, 1), (-1, 0, 0)])
#    degenerate_equal(
#        crossmap.coding_to_coordinate, [(10, 1, 1), (11, 0, 1), (1, 0, 2)])
#
#
#def test_Crossmap_inverted_no_utr_degenerate():
#    crossmap = Crossmap([(10, 20)], (10, 20), True)
#
#    degenerate_equal(
#        crossmap.coding_to_coordinate, [(1, -1, 1), (-1, 0, 1), (-1, 0, 0)])
#    degenerate_equal(
#        crossmap.coding_to_coordinate, [(10, 1, 1), (11, 0, 1), (1, 0, 2)])
#
#
#def test_Crossmap_degenerate_return():
#    crossmap = Crossmap([(10, 20)], (11, 19))
#
#    assert crossmap.coordinate_to_coding(9, True) == (-2, 0, 0)
#    assert crossmap.coordinate_to_coding(20, True) == (2, 0, 2)
#
#
#def test_Crossmap_inverted_degenerate_return():
#    crossmap = Crossmap([(10, 20)], (11, 19), True)
#
#    assert crossmap.coordinate_to_coding(9, True) == (2, 0, 2)
#    assert crossmap.coordinate_to_coding(20, True) == (-2, 0, 0)
#
#
#def test_Crossmap_no_utr_degenerate_return():
#    crossmap = Crossmap([(10, 20), (30, 40)], (10, 40))
#
#    assert crossmap.coordinate_to_coding(9, True) == (-1, 0, 0)
#    assert crossmap.coordinate_to_coding(40, True) == (1, 0, 2)
#
#
#def test_Crossmap_no_utr_inverted_degenerate_return():
#    crossmap = Crossmap([(10, 20), (30, 40)], (10, 40), True)
#
#    assert crossmap.coordinate_to_coding(9, True) == (1, 0, 2)
#    assert crossmap.coordinate_to_coding(40, True) == (-1, 0, 0)
#
#
#def test_Crossmap_degenerate_cds_return():
#    crossmap = Crossmap([(10, 11)], (10, 11))
#
#    assert crossmap.coordinate_to_coding(9, True) == (-1, 0, 0)
#    assert crossmap.coordinate_to_coding(11, True) == (1, 0, 2)
#
#
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
