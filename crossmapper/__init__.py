"""
crossmapper: Mutalyzer crossmapper.


Copyright (c) 2016 Jeroen F.J. Laros <J.F.J.Laros@lumc.nl>
Copyright (c) 2016 Leiden University Medical Center <humgen@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""
__version_info__ = ('0', '0', '1')

__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Jeroen F.J. Laros'
__contact__ = 'J.F.J.Laros@lumc.nl'
__homepage__ = 'https://github.com/mutalyzer/crossmapper'

usage = __doc__.split('\n\n\n')

from . import crossmapper

Crossmap = crossmapper.Crossmap

def doc_split(func):
    return func.__doc__.split('\n\n')[0]


def version(name):
    return '%s version %s\n\nAuthor   : %s <%s>\nHomepage : %s' % (name,
        __version__, __author__, __contact__, __homepage__)
