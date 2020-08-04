from mutalyzer_crossmapper import Locus

from helper import degenerate_equal, invariant


def test_Locus():
    locus = Locus((30, 35))

    invariant(locus.to_position, 29, locus.to_coordinate, (1, -1))
    invariant(locus.to_position, 30, locus.to_coordinate, (1, 0))
    invariant(locus.to_position, 31, locus.to_coordinate, (2, 0))
    invariant(locus.to_position, 33, locus.to_coordinate, (4, 0))
    invariant(locus.to_position, 34, locus.to_coordinate, (5, 0))
    invariant(locus.to_position, 35, locus.to_coordinate, (5, 1))


def test_Locus_inverted():
    locus = Locus((30, 35), True)

    invariant(locus.to_position, 35, locus.to_coordinate, (1, -1))
    invariant(locus.to_position, 34, locus.to_coordinate, (1, 0))
    invariant(locus.to_position, 33, locus.to_coordinate, (2, 0))
    invariant(locus.to_position, 31, locus.to_coordinate, (4, 0))
    invariant(locus.to_position, 30, locus.to_coordinate, (5, 0))
    invariant(locus.to_position, 29, locus.to_coordinate, (5, 1))


def test_Locus_degenerate():
    locus = Locus((10, 20))

    degenerate_equal(locus.to_coordinate, 9, [(1, -1), (-1, 0)])
    degenerate_equal(locus.to_coordinate, 20, [(10, 1), (11, 0)])


def test_Locus_inverted_degenerate():
    locus = Locus((10, 20), True)

    degenerate_equal(locus.to_coordinate, 20, [(1, -1), (-1, 0)])
    degenerate_equal(locus.to_coordinate, 9, [(10, 1), (11, 0)])
