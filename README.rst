HGVS position crossmapper
=========================

.. image:: https://img.shields.io/github/last-commit/mutalyzer/crossmapper.svg
   :target: https://github.com/mutalyzer/crossmapper/graphs/commit-activity
.. image:: https://github.com/mutalyzer/crossmapper/actions/workflows/python-package.yml/badge.svg
   :target: https://github.com/mutalyzer/crossmapper/actions/workflows/python-package.yml
.. image:: https://readthedocs.org/projects/mutalyzer-crossmapper/badge/?version=latest
   :target: https://mutalyzer-crossmapper.readthedocs.io/en/latest
.. image:: https://img.shields.io/github/release-date/mutalyzer/crossmapper.svg
   :target: https://github.com/mutalyzer/crossmapper/releases
.. image:: https://img.shields.io/github/release/mutalyzer/crossmapper.svg
   :target: https://github.com/mutalyzer/crossmapper/releases
.. image:: https://img.shields.io/pypi/v/mutalyzer-crossmapper.svg
   :target: https://pypi.org/project/mutalyzer-crossmapper/
.. image:: https://img.shields.io/github/languages/code-size/mutalyzer/crossmapper.svg
   :target: https://github.com/mutalyzer/crossmapper
.. image:: https://img.shields.io/github/languages/count/mutalyzer/crossmapper.svg
   :target: https://github.com/mutalyzer/crossmapper
.. image:: https://img.shields.io/github/languages/top/mutalyzer/crossmapper.svg
   :target: https://github.com/mutalyzer/crossmapper
.. image:: https://img.shields.io/github/license/mutalyzer/crossmapper.svg
   :target: https://raw.githubusercontent.com/mutalyzer/crossmapper/master/LICENSE.md

----

This library provides an interface to convert (cross map) between different
HGVS numbering_ systems.

Converting between the transcript oriented c. or n. and the genomic oriented g.
numbering systems can be difficult, especially when the transcript in question
resides on the complement strand.

**Features:**

- Support for genomic positions to standard coordinates and vice versa.
- Support for noncoding positions to standard coordinates and vice versa.
- Support for coding positions to standard coordinates and vice versa.
- Support for protein positions to standard coordinates and vice versa.
- Basic classes for loci that can be used for genomic loci other than genes.

Please see ReadTheDocs_ for the latest documentation.

Quick start
-----------
The `Genomic` class provides an interface to conversions between genomic positions and coordinates.

**Genomic Position Model**
Genomic positions follow the HGVS ``g`` coordinate system. They are represented
as dictionaries:

.. code:: json

    {"position": int}

Where:

- **position**: a positive integer

**Genomic Position Conversion**
.. code:: python

    >>> from mutalyzer_crossmapper import Genomic
    >>> crossmap = Genomic()
    >>> crossmap.coordinate_to_genomic(0)
    {"position": 1}
    >>> crossmap.genomic_to_coordinate({"position": 1})
    0

On top of the functionality provided by the ``Genomic`` class, the
``NonCoding`` class provides an interface to conversions between noncoding
positions and coordinates.

**NonCoding Position Model**
Noncoding positions follow the HGVS `n` coordinate system. They are represented as dictionaries:
.. code:: json
   {
    "position": int,
    "offset": int,
    "region": str
   }
Where:

- **position**: a positive interger
- **offset**: an interger indicating the offset relative to the position (e.g.,
negative for upstream or positive for downstream)
- **region**: a string describing the region type (empty string `""` for standard
noncoding positions, `"u"` for upstream and `"d"` for downstream.)

**NonCoding Position Conversion**
.. code:: python

    >>> from mutalyzer_crossmapper import NonCoding
    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> crossmap = NonCoding(exons)
    >>> crossmap.coordinate_to_noncoding(35)
    {"position": 14, "offset": 1, "region": ""}
    >>> crossmap.noncoding_to_coordinate({"position": 14, "offset": 1, "region": ""})
    35

**Notes**
- Add the flag ``inverted=True`` to the constructor when the transcript resides
on the reverse complement strand.

On top of the functionality provided by the ``NonCoding`` class, the ``Coding``
class provides an interface to conversions between coding positions and
coordinates as well as conversions between protein positions and coordinates.

**Coding Position Model**
Coding positions follow the HGVS `c`` coordinate system. They are represented as
dictionaries:
.. code:: json
   {
    "position": int,
    "offset": int,
    "region": str
   }
Where:

- **position**: a positive interger
- **offset**: an interger indicating the offset relative to the position (e.g.,
negative for upstream or positive for downstream)
- **region**: a string describing the region type (empty string `""` for standard
coding positions, `"-"` for 5' UTR, `"*"` for 3' UTR, `"u"` for upstream and `"d"`
for downstream.)

**Coding Position Conversion**
.. code:: python

    >>> from mutalyzer_crossmapper import Coding
    >>> cds = (32, 43)
    >>> crossmap = Coding(exons, cds)
    >>> crossmap.coordinate_to_coding(31)
    {"position": -1, "offset": 0, "region": "-"}
    >>> crossmap.coding_to_coordinate({"position": -1, "offset": 0, "region": "-"})
    31

**Notes**
- Again, the flag ``inverted=True`` can be used for transcripts that reside on the reverse complement strand.

**Protein Position Model**
Protein positions follow the HGVS `p`` coordinate system. They are represented
as dictionaries:
.. code:: json
{
    "position": int,
    "position_in_codon": int,
    "offset": int,
    "region": str
}
Where:

- **position**: the amino acid position (1-based)
- **position_in_codon**: the codon nucleotide index (1, 2, or 3)
- **offset**: an integer indicating offset relative to the codon
- **region**: a string describing the region type (empty string `""`` for standard positions)

**Protein Position Conversion**

Conversions between protein positions and coordinates are done as follows.
.. code:: python

    >>> crossmap.coordinate_to_protein(41)
    {"position": 2, "position_in_codon": 2, "offset": 1, "region": ""}
    >>> crossmap.protein_to_coordinate({"position": 2, "position_in_codon": 2, "offset": 1, "region": ""})
    41


.. _numbering: http://varnomen.hgvs.org/bg-material/numbering/
.. _ReadTheDocs: https://mutalyzer-crossmapper.readthedocs.io
