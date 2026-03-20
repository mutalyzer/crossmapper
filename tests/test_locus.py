from mutalyzer_crossmapper import Locus

from helper import degenerate_equal, invariant


def test_Locus():
    """Forward orientent Lovus."""
    locus = Locus((30, 35))

    invariant(locus.to_position, 29, locus.to_coordinate, {'position': 0, 'offset': -1})
    invariant(locus.to_position, 30, locus.to_coordinate, {'position': 0, 'offset': 0})
    invariant(locus.to_position, 31, locus.to_coordinate, {'position': 1, 'offset': 0})
    invariant(locus.to_position, 33, locus.to_coordinate, {'position': 3, 'offset': 0})
    invariant(locus.to_position, 34, locus.to_coordinate, {'position': 4, 'offset': 0})
    invariant(locus.to_position, 35, locus.to_coordinate, {'position': 4, 'offset': 1})


def test_Locus_inverted():
    """Reverse orientent Lovus."""
    locus = Locus((30, 35), True)

    invariant(locus.to_position, 35, locus.to_coordinate, {'position': 0, 'offset': -1})
    invariant(locus.to_position, 34, locus.to_coordinate, {'position': 0, 'offset': 0})
    invariant(locus.to_position, 33, locus.to_coordinate, {'position': 1, 'offset': 0})
    invariant(locus.to_position, 31, locus.to_coordinate, {'position': 3, 'offset': 0})
    invariant(locus.to_position, 30, locus.to_coordinate, {'position': 4, 'offset': 0})
    invariant(locus.to_position, 29, locus.to_coordinate, {'position': 4, 'offset': 1})


def test_Locus_degenerate():
    """Degenerate positions are silently corrected."""
    locus = Locus((10, 20))

    degenerate_equal(locus.to_coordinate, 9, [{'position': 0, 'offset': -1}, {'position': -1, 'offset': 0}])
    degenerate_equal(locus.to_coordinate, 20, [{'position': 9, 'offset': 1}, {'position': 10, 'offset': 0}])


def test_Locus_inverted_degenerate():
    """Degenerate positions are silently corrected."""
    locus = Locus((10, 20), True)

    degenerate_equal(locus.to_coordinate, 20, [{'position': 0, 'offset': -1}, {'position': -1, 'offset': 0}])
    degenerate_equal(locus.to_coordinate, 9, [{'position': 9, 'offset': 1}, {'position': 10, 'offset': 0}])
