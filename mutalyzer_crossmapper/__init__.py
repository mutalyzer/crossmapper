"""Crossmapper position conversion library.

Definitions:

- Coordinates are zero based, non-negative integers.
- Locations are zero based right-open non-negative integer intervals,
  consistent with Python's range() and sequence slicing functions.
- Loci and exons are locations.
- An exon list is a list of locations that, when flattened, is an increasing
  sequence.
- A position is a 2-tuple of which the first element is a one based non-zero
  integer relative to an element in a location and the second element is an
  integer offset relative to the first element.
"""
from pkg_resources import get_distribution

from .crossmapper import Coding, Genomic, NonCoding
from .location import nearest_location
from .locus import Locus
from .multi_locus import MultiLocus


def _get_metadata(name):
    pkg = get_distribution(__package__)

    for line in pkg.get_metadata_lines(pkg.PKG_INFO):
        if line.startswith('{}: '.format(name)):
            return line.split(': ')[1]

    return ''


_copyright_notice = 'Copyright (c) {} <{}>'.format(
    _get_metadata('Author'), _get_metadata('Author-email'))

usage = [_get_metadata('Summary'), _copyright_notice]


def doc_split(func):
    return func.__doc__.split('\n\n')[0]


def version(name):
    return '{} version {}\n\n{}\nHomepage: {}'.format(
        _get_metadata('Name'), _get_metadata('Version'), _copyright_notice,
        _get_metadata('Home-page'))
