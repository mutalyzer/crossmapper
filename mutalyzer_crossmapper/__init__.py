from importlib.metadata import metadata

from .crossmapper import Coding, Genomic, NonCoding
from .location import nearest_location
from .locus import Locus
from .multi_locus import MultiLocus
from .models import Point, GenomicPoint, NonCodingPoint, CodingPoint, ProteinPoint


def _get_metadata(name: str) -> str:
    """Get metadata from the package using importlib.metadata"""
    try:
        meta = metadata('mutalyzer_crossmapper')
        return meta.get(name, '')
    except Exception:
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
