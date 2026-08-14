from mutalyzer_crossmapper.location import nearest_location


def _check_int(value: int) -> None:
    """Check if the value is a non-negative integer.

    :arg int value: Value to check.

    :raises ValueError: If the value is invalid.
    """
    if not isinstance(value, int):
        raise ValueError("Value must be an integer.")


def _check_in_range(value: int, length: int) -> None:
    if value > length:
        raise ValueError(f"Value {value} must be within the bounds of the reference sequence {length}.")


def _check_non_negative(value: int, length: int|None = None) -> None:
    """Check if the coordinate is a non-negative integer.

    :arg int value: Value to check.

    :raises ValueError: If the coordinate is invalid.
    """
    _check_int(value)
    if value < 0:
        raise ValueError("Value must be non-negative.")
    if length is not None:
        _check_in_range(value, length)


def _check_locus(locus: tuple[int, int], length: int| None = None) -> None:
    """Check if the range is valid.

    :arg tuple[int, int] locus: Locus to check.

    :raises ValueError: If the range is invalid.
    """
    if len(locus) != 2:
        raise ValueError("Locus must be a tuple of two values.")

    for value in locus:
        _check_non_negative(value, length)

    if locus[0] > locus[1]:
        raise ValueError("Start of locus must be smaller than or equal to end of locus.")


def _check_exons(exons: list[tuple[int, int]], length: int|None = None) -> None:
    """Check if the exons are valid.
    The exons are valid
    if they are a list of valid loci,
    non-overlapping,
    and within the bounds of the reference sequence.

    :arg list[tuple[int, int]] exons: Exons to check.

    :raises ValueError: If the exons are invalid.
    """
    for exon in exons:
        _check_locus(exon, length)

    for e1, e2 in zip(exons, exons[1:]):
        if e2[0] < e1[1]:
            raise ValueError(f"Exon {e2} and exon {e1} are overlapping.")


def _check_cds(cds: tuple[int, int], exons: list[tuple[int, int]], length: int|None = None) -> None:
    """Check if the CDS is valid.

    :arg tuple[int, int] cds: CDS to check.
    :arg list[tuple[int, int]] exons: List of exons.
    :arg int|None length: Length of the reference sequence.

    :raises ValueError: If the CDS is invalid.
    """
    _check_locus(cds, length)
    for coord in cds:
        index = nearest_location(exons, coord)
        if coord < exons[index][0] or coord >= exons[index][1]:
            raise ValueError(f"Coordinate {coord} of CDS {cds} is not within any exon.")

