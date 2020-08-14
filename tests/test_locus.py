from mutalyzer_crossmapper import Locus

from helper import degenerate_equal, invariant


def test_Locus():
    """Forward orientent Lovus."""
    locus = Locus((30, 35))

    invariant(locus.to_position, 29, locus.to_coordinate, (0, -1))
    invariant(locus.to_position, 30, locus.to_coordinate, (0, 0))
    invariant(locus.to_position, 31, locus.to_coordinate, (1, 0))
    invariant(locus.to_position, 33, locus.to_coordinate, (3, 0))
    invariant(locus.to_position, 34, locus.to_coordinate, (4, 0))
    invariant(locus.to_position, 35, locus.to_coordinate, (4, 1))


def test_Locus_inverted():
    """Reverse orientent Lovus."""
    locus = Locus((30, 35), True)

    invariant(locus.to_position, 35, locus.to_coordinate, (0, -1))
    invariant(locus.to_position, 34, locus.to_coordinate, (0, 0))
    invariant(locus.to_position, 33, locus.to_coordinate, (1, 0))
    invariant(locus.to_position, 31, locus.to_coordinate, (3, 0))
    invariant(locus.to_position, 30, locus.to_coordinate, (4, 0))
    invariant(locus.to_position, 29, locus.to_coordinate, (4, 1))


def test_Locus_degenerate():
    """Degenerate positions are silently corrected."""
    locus = Locus((10, 20))

    degenerate_equal(locus.to_coordinate, 9, [(0, -1), (-1, 0)])
    degenerate_equal(locus.to_coordinate, 20, [(9, 1), (10, 0)])


#def test_Locus_degenerate_return():
#    """Degenerate position can be retured."""
#    locus = Locus((10, 20))
#
#    assert locus.to_position(9, True) == (-1, 0)
#    assert locus.to_position(20, True) == (10, 0)


def test_Locus_inverted_degenerate():
    """Degenerate positions are silently corrected."""
    locus = Locus((10, 20), True)

    degenerate_equal(locus.to_coordinate, 20, [(0, -1), (-1, 0)])
    degenerate_equal(locus.to_coordinate, 9, [(9, 1), (10, 0)])


#def test_Locus_inverted_degenerate_return():
#    """Degenerate position can be retured."""
#    locus = Locus((10, 20), True)
#
#    assert locus.to_position(20, True) == (-1, 0)
#    assert locus.to_position(9, True) == (11, 0)
