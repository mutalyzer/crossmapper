from mutalyzer_crossmapper import Coding, Genomic, NonCoding

from helper import degenerate_equal, invariant

_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_cds = (32, 43)


def test_Genomic():
    """Genomic positions are coordinates incremented by one."""
    crossmap = Genomic()

    invariant(
        crossmap.coordinate_to_genomic,
        0,
        crossmap.genomic_to_coordinate,
        {'position': 1},
    )
    invariant(
        crossmap.coordinate_to_genomic,
        98,
        crossmap.genomic_to_coordinate,
        {'position': 99},
    )


def test_NonCoding():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons)

    # Boundary between upstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding,
        3,
        crossmap.noncoding_to_coordinate,
        {'position': 2, 'offset': 0, 'region': 'u'},
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        4,
        crossmap.noncoding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'u'},
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        5,
        crossmap.noncoding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )

    # Boundary between downstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding,
        71,
        crossmap.noncoding_to_coordinate,
        {'position': 22, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        72,
        crossmap.noncoding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'd'},
    )


def test_NonCoding_inverted():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons, inverted=True)

    # Boundary between upstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding,
        72,
        crossmap.noncoding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'u'},
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        71,
        crossmap.noncoding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )

    # Boundary between downstream and transcript.
    invariant(
        crossmap.coordinate_to_noncoding,
        5,
        crossmap.noncoding_to_coordinate,
        {'position': 22, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_noncoding,
        4,
        crossmap.noncoding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'd'},
    )


def test_NonCoding_degenerate():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons)

    # Boundary between upstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate,
        4,
        [
            {'position': 1, 'offset': -1, 'region': ''},
            {'position': 1, 'offset': 0, 'region': 'u'},
            {'position': 1, 'offset': 0, 'region': '-'},
            {'position': 2, 'offset': 1, 'region': '-'},
        ],
    )

    # Boundary between downstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate,
        72,
        [
            {'position': 1, 'offset': 0, 'region': 'd'},
            {'position': 22, 'offset': 1, 'region': ''},
            {'position': 23, 'offset': 0, 'region': ''},
            {'position': 24, 'offset': -1, 'region': ''},
            {'position': 1, 'offset': 0, 'region': '*'},
        ],
    )


def test_NonCoding_inverted_degenerate():
    """Forward oriented noncoding transcript."""
    crossmap = NonCoding(_exons, inverted=True)

    # Boundary between upstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate,
        72,
        [
            {'position': 1, 'offset': -1, 'region': ''},
            {'position': 1, 'offset': 0, 'region': 'u'},
            {'position': 1, 'offset': 0, 'region': '-'},
        ],
    )

    # Boundary between downstream and transcript.
    degenerate_equal(
        crossmap.noncoding_to_coordinate,
        4,
        [
            {'position': 1, 'offset': 0, 'region': 'd'},
            {'position': 1, 'offset': 0, 'region': '*'},
            {'position': 23, 'offset': 0, 'region': ''},
            {'position': 22, 'offset': 1, 'region': ''},
        ],
    )


def test_NonCoding_degenerate_return():
    crossmap = NonCoding(_exons)

    assert crossmap.coordinate_to_noncoding(4, True) == {
        'position': 1,
        'offset': 0,
        'region': '-',
    }

    assert crossmap.coordinate_to_noncoding(72, True) == {
        'position': 1,
        'offset': 0,
        'region': '*',
    }


def test_NonCoding_inverted_degenerate_return():
    crossmap = NonCoding(_exons, True)

    assert crossmap.coordinate_to_noncoding(72, True) == {
        'position': 1,
        'offset': 0,
        'region': '-',
    }

    assert crossmap.coordinate_to_noncoding(4, True) == {
        'position': 1,
        'offset': 0,
        'region': '*',
    }


def test_NonCoding_degenerate_no_return():
    """Degenerate internal positions do not exist."""
    crossmap = NonCoding(_exons)

    assert crossmap.coordinate_to_noncoding(25) == crossmap.coordinate_to_noncoding(25, True)


def test_NonCoding_inverted_degenerate_no_return():
    """Degenerate internal positions do not exist."""
    crossmap = NonCoding(_exons, True)

    assert crossmap.coordinate_to_noncoding(25) == crossmap.coordinate_to_noncoding(25, True)


def test_Coding():
    """Forward oriented coding transcript."""
    crossmap = Coding(_exons, _cds)

    # Boundary between 5' and CDS.
    invariant(
        crossmap.coordinate_to_coding,
        31,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': '-'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        32,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )

    # Boundary between CDS and 3'.
    invariant(
        crossmap.coordinate_to_coding,
        42,
        crossmap.coding_to_coordinate,
        {'position': 6, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_coding,
        43,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': '*'},
    )


def test_Coding_inverted():
    """Reverse oriented coding transcript."""
    crossmap = Coding(_exons, _cds, True)

    # Boundary between 5' and CDS.
    invariant(
        crossmap.coordinate_to_coding,
        43,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': '-'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        42,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )

    # Boundary between CDS and 3'.
    invariant(
        crossmap.coordinate_to_coding,
        32,
        crossmap.coding_to_coordinate,
        {'position': 6, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_coding,
        31,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': '*'},
    )


def test_Coding_regions():
    """The CDS can start or end on a region boundary."""
    crossmap = Coding([(10, 21), (30, 40), (49, 60)], (30, 40))

    # Upstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding,
        25,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 5, 'region': '-'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        26,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': -4, 'region': ''},
    )

    # Downstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding,
        44,
        crossmap.coding_to_coordinate,
        {'position': 10, 'offset': 5, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_coding,
        45,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': -4, 'region': '*'},
    )


def test_Coding_regions_inverted():
    """The CDS can start or end on a region boundary."""
    crossmap = Coding([(10, 21), (30, 40), (49, 60)], (30, 40), True)

    # Upstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding,
        44,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 5, 'region': '-'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        43,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': -4, 'region': ''},
    )

    # Downstream odd length intron between two regions.
    invariant(
        crossmap.coordinate_to_coding,
        25,
        crossmap.coding_to_coordinate,
        {'position': 10, 'offset': 5, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_coding,
        24,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': -4, 'region': '*'},
    )


def test_Coding_no_utr5():
    """A 5' UTR may be missing."""
    crossmap = Coding([(10, 20)], (10, 15))

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding,
        9,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'u'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        10,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )


def test_Coding_no_utr5_inverted():
    """A 5' UTR may be missing."""
    crossmap = Coding([(10, 20)], (15, 20), True)

    # Direct transition from upstream to CDS.
    invariant(
        crossmap.coordinate_to_coding,
        20,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'u'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        19,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )


def test_Coding_no_utr3():
    """A 3' UTR may be missing."""
    crossmap = Coding([(10, 20)], (15, 20))

    # Direct transition from CDS to downstream.
    invariant(
        crossmap.coordinate_to_coding,
        19,
        crossmap.coding_to_coordinate,
        {'position': 5, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_coding,
        20,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'd'},
    )


def test_Coding_no_utr3_inverted():
    """A 3' UTR may be missing."""
    crossmap = Coding([(10, 20)], (10, 15), True)

    # Direct transition from CDS to downstream.
    invariant(
        crossmap.coordinate_to_coding,
        10,
        crossmap.coding_to_coordinate,
        {'position': 5, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_coding,
        9,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'd'},
    )


def test_Coding_small_utr5():
    """A 5' UTR may be of lenght one."""
    crossmap = Coding([(10, 20)], (11, 15))

    # Transition from upstream to 5' UTR to CDS.
    invariant(
        crossmap.coordinate_to_coding,
        9,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'u'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        10,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': '-'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        11,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )


def test_Coding_small_utr5_inverted():
    """A 5' UTR may be of lenght one."""
    crossmap = Coding([(10, 20)], (15, 19), True)

    # Transition from upstream to 5' UTR to CDS.
    invariant(
        crossmap.coordinate_to_coding,
        20,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'u'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        19,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': '-'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        18,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': ''},
    )


def test_Coding_small_utr3():
    """A 5' UTR may be of lenght one."""
    crossmap = Coding([(10, 20)], (15, 19))

    # Transition from CDS to 3' UTR to downstream.
    invariant(
        crossmap.coordinate_to_coding,
        18,
        crossmap.coding_to_coordinate,
        {'position': 4, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_coding,
        19,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': '*'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        20,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'd'},
    )


def test_Coding_small_utr3_inverted():
    """A 5' UTR may be of lenght one."""
    crossmap = Coding([(10, 20)], (11, 15), True)

    # Transition from CDS to 3' UTR to downstream.
    invariant(
        crossmap.coordinate_to_coding,
        11,
        crossmap.coding_to_coordinate,
        {'position': 4, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_coding,
        10,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': '*'},
    )
    invariant(
        crossmap.coordinate_to_coding,
        9,
        crossmap.coding_to_coordinate,
        {'position': 1, 'offset': 0, 'region': 'd'},
    )


def test_Coding_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    crossmap = Coding([(10, 20)], (11, 19))

    # Degenerate position in upstream.
    degenerate_equal(
        crossmap.coding_to_coordinate,
        9,
        [
            {'position': 1, 'offset': 0, 'region': 'u'},
            {'position': 2, 'offset': 0, 'region': '-'},
            {'position': 1, 'offset': -2, 'region': ''},
            {'position': 1, 'offset': -10, 'region': '*'},
            {'position': 2, 'offset': -11, 'region': '*'},
            {'position': 3, 'offset': 1, 'region': '-'},
            {'position': 4, 'offset': 2, 'region': '-'},
        ],
    )
    degenerate_equal(
        crossmap.coding_to_coordinate,
        20,
        [
            {'position': 1, 'offset': 0, 'region': 'd'},
            {'position': 2, 'offset': 0, 'region': '*'},
            {'position': 8, 'offset': 2, 'region': ''},
            {'position': 1, 'offset': 10, 'region': '-'},
            {'position': 2, 'offset': 11, 'region': '-'},
            {'position': 7, 'offset': 3, 'region': ''},
        ],
    )


def test_Coding_inverted_degenerate():
    """Degenerate upstream and downstream positions are silently corrected."""
    crossmap = Coding([(10, 20)], (11, 19), True)

    degenerate_equal(
        crossmap.coding_to_coordinate,
        20,
        [
            {'position': 1, 'offset': 0, 'region': 'u'},
            {'position': 2, 'offset': 0, 'region': '-'},
            {'position': 1, 'offset': -2, 'region': ''},
            {'position': 1, 'offset': -10, 'region': '*'},
            {'position': 1, 'offset': -11, 'region': 'd'},
            {'position': 2, 'offset': -3, 'region': ''},
        ],
    )
    degenerate_equal(
        crossmap.coding_to_coordinate,
        9,
        [
            {'position': 1, 'offset': 0, 'region': 'd'},
            {'position': 2, 'offset': 0, 'region': '*'},
            {'position': 8, 'offset': 2, 'region': ''},
            {'position': 1, 'offset': 10, 'region': '-'},
            {'position': 1, 'offset': 11, 'region': 'u'},
            {'position': 2, 'offset': 12, 'region': 'u'},

        ],
    )


def test_Coding_degenerate_return():
    """Degenerate upstream and downstream positions may be returned."""
    crossmap = Coding([(10, 20)], (11, 19))

    assert crossmap.coordinate_to_coding(9, True) == {
        'position': 2,
        'offset': 0,
        'region': '-',
    }
    assert crossmap.coordinate_to_coding(20, True) == {
        'position': 2,
        'offset': 0,
        'region': '*',
    }


def test_Coding_inverted_degenerate_return():
    """Degenerate upstream and downstream positions may be returned."""
    crossmap = Coding([(10, 20)], (11, 19), True)

    assert crossmap.coordinate_to_coding(20, True) == {
        'position': 2,
        'offset': 0,
        'region': '-',
    }
    assert crossmap.coordinate_to_coding(9, True) == {
        'position': 2,
        'offset': 0,
        'region': '*',
    }


def test_Coding_degenerate_no_return():
    """Degenerate internal positions do not exist."""
    crossmap = Coding([(10, 20), (30, 40)], (10, 40))

    assert crossmap.coordinate_to_coding(25) == crossmap.coordinate_to_coding(25, True)


def test_Coding_inverted_degenerate_no_return():
    """Degenerate internal positions do not exist."""
    crossmap = Coding([(10, 20), (30, 40)], (10, 40), True)

    assert crossmap.coordinate_to_coding(25) == crossmap.coordinate_to_coding(25, True)


def test_Coding_no_utr_degenerate():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11))

    degenerate_equal(
        crossmap.coding_to_coordinate,
        9,
        [
            {'position': 1, 'offset': 0, 'region': 'u'},
            {'position': 1, 'offset': 0, 'region': '-'},
            {'position': 1, 'offset': -2, 'region': '*'},
            {'position': 1, 'offset': -1, 'region': ''},
            {'position': 1, 'offset': -2, 'region': 'd'},
        ],
    )
    degenerate_equal(
        crossmap.coding_to_coordinate,
        11,
        [
            {'position': 1, 'offset': 0, 'region': 'd'},
            {'position': 1, 'offset': 0, 'region': '*'},
            {'position': 1, 'offset': 2, 'region': '-'},
            {'position': 1, 'offset': 1, 'region': ''},
            {'position': 1, 'offset': 2, 'region': 'u'},
        ],
    )


def test_Coding_inverted_no_utr_degenerate():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11), True)

    degenerate_equal(
        crossmap.coding_to_coordinate,
        11,
        [
            {'position': 1, 'offset': 0, 'region': 'u'},
            {'position': 1, 'offset': 0, 'region': '-'},
            {'position': 1, 'offset': -2, 'region': '*'},
            {'position': 1, 'offset': -1, 'region': ''},
            {'position': 1, 'offset': -2, 'region': 'd'},
        ],
    )
    degenerate_equal(
        crossmap.coding_to_coordinate,
        9,
        [
            {'position': 1, 'offset': 0, 'region': 'd'},
            {'position': 1, 'offset': 0, 'region': '*'},
            {'position': 1, 'offset': 2, 'region': '-'},
            {'position': 1, 'offset': 1, 'region': ''},
            {'position': 1, 'offset': 2, 'region': 'u'},
        ],
    )


def test_Coding_no_utr_degenerate_return():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11))

    assert crossmap.coordinate_to_coding(8, True) == {
        'position': 2,
        'offset': 0,
        'region': '-',
    }
    assert crossmap.coordinate_to_coding(9, True) == {
        'position': 1,
        'offset': 0,
        'region': '-',
    }
    assert crossmap.coordinate_to_coding(11, True) == {
        'position': 1,
        'offset': 0,
        'region': '*',
    }
    assert crossmap.coordinate_to_coding(12, True) == {
        'position': 2,
        'offset': 0,
        'region': '*',
    }


def test_Coding_inverted_no_utr_degenerate_return():
    """UTRs may be missing."""
    crossmap = Coding([(10, 11)], (10, 11), True)

    assert crossmap.coordinate_to_coding(11, True) == {
        'position': 1,
        'offset': 0,
        'region': '-',
    }
    assert crossmap.coordinate_to_coding(9, True) == {
        'position': 1,
        'offset': 0,
        'region': '*',
    }


def test_Coding_protein():
    """Protein positions."""
    crossmap = Coding(_exons, _cds)

    # Boundary between upstream and 5' UTR
    invariant(
        crossmap.coordinate_to_protein,
        4,
        crossmap.protein_to_coordinate,
        {'position': 1, 'position_in_codon': 3, 'offset': 0, 'region': 'u'}
    )
    invariant(
        crossmap.coordinate_to_protein,
        5,
        crossmap.protein_to_coordinate,
        {'position': 4, 'position_in_codon': 2, 'offset': 0, 'region': '-'}
    )

    # Boundary between 5' UTR and CDS
    invariant(
        crossmap.coordinate_to_protein,
        31,
        crossmap.protein_to_coordinate,
        {'position': 1, 'position_in_codon': 3, 'offset': 0, 'region': '-'},
    )
    invariant(
        crossmap.coordinate_to_protein,
        32,
        crossmap.protein_to_coordinate,
        {'position': 1, 'position_in_codon': 1, 'offset': 0, 'region': ''},
    )

    # Intron boundary.
    invariant(
        crossmap.coordinate_to_protein,
        34,
        crossmap.protein_to_coordinate,
        {'position': 1, 'position_in_codon': 3, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_protein,
        35,
        crossmap.protein_to_coordinate,
        {'position': 1, 'position_in_codon': 3, 'offset': 1, 'region': ''},
    )

    # Boundary between CDS and 3' UTR.
    invariant(
        crossmap.coordinate_to_protein,
        42,
        crossmap.protein_to_coordinate,
        {'position': 2, 'position_in_codon': 3, 'offset': 0, 'region': ''},
    )
    invariant(
        crossmap.coordinate_to_protein,
        43,
        crossmap.protein_to_coordinate,
        {'position': 1, 'position_in_codon': 1, 'offset': 0, 'region': '*'},
    )

    # Boundary between 3' UTR and downstream
    invariant(
        crossmap.coordinate_to_protein,
        71,
        crossmap.protein_to_coordinate,
        {'position': 2, 'position_in_codon': 2, 'offset': 0, 'region': '*'}
    )
    invariant(
        crossmap.coordinate_to_protein,
        72,
        crossmap.protein_to_coordinate,
        {'position': 1, 'position_in_codon': 1, 'offset': 0, 'region': 'd'}
    )
