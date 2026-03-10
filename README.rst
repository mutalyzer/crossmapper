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

Quick Start
===========

An example below uses the following transcript data:

.. code-block:: python

    >>>_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>>_cds = (32, 43)


Genomic Class
-------------

The ``Genomic`` class provides an interface for conversions between genomic positions and coordinates.

Genomic Position Model
~~~~~~~~~~~~~~~~~~~~~~~

Genomic positions follow the HGVS ``g`` coordinate system. They are represented as dictionaries:

.. code-block:: json

    {
        "position": int
    }

Where:

- **position**: a positive integer

Genomic Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from mutalyzer_crossmapper import Genomic
    >>> crossmap = Genomic()
    >>> crossmap.coordinate_to_genomic(0)
    {"position": 1}
    >>> crossmap.genomic_to_coordinate({"position": 1})
    0

Here is the mapping of coordinates to genomic positions:
+------------+----------+
| Coordinate | Position |
+============+==========+
| 0          | 1        |
+------------+----------+
| 1          | 2        |
+------------+----------+
| 2          | 3        |
+------------+----------+
| 3          | 4        |
+------------+----------+
| 4          | 5        |
+------------+----------+
| 5          | 6        |
+------------+----------+
| 6          | 7        |
+------------+----------+
| 7          | 8        |
+------------+----------+
| 8          | 9        |
+------------+----------+
| 9          | 10       |
+------------+----------+
| 10         | 11       |
+------------+----------+
| 11         | 12       |
+------------+----------+
| 12         | 13       |
+------------+----------+
| 13         | 14       |
+------------+----------+
| 14         | 15       |
+------------+----------+
| 15         | 16       |
+------------+----------+
| 16         | 17       |
+------------+----------+
| 17         | 18       |
+------------+----------+
| 18         | 19       |
+------------+----------+
| 19         | 20       |
+------------+----------+

NonCoding Class
---------------

The ``NonCoding`` class provides conversions between noncoding positions and coordinates.

NonCoding Position Model
~~~~~~~~~~~~~~~~~~~~~~~

Noncoding positions follow the HGVS ``n`` coordinate system. They are represented as dictionaries:

.. code-block:: json

    {
        "position": int,
        "offset": int,
        "region": str
    }

Where:

- **position**: a positive integer
- **offset**: an integer indicating the offset relative to the position (negative for upstream, positive for downstream)
- **region**: a string describing the region type (``""`` for standard, ``"u"`` for upstream, ``"d"`` for downstream)

NonCoding Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from mutalyzer_crossmapper import NonCoding
    >>> crossmap = NonCoding(_exons)
    >>> crossmap.coordinate_to_noncoding(35)
    {"position": 14, "offset": 1, "region": ""}
    >>> crossmap.noncoding_to_coordinate({"position": 14, "offset": 1, "region": ""})
    35

Notes
~~~~~

- Add the flag ``inverted=True`` to the constructor when the transcript resides on the reverse complement strand.

Coding Class
------------

The ``Coding`` class provides conversions between coding positions and coordinates, as well as protein positions.

Coding Position Model
~~~~~~~~~~~~~~~~~~~~

Coding positions follow the HGVS ``c`` coordinate system. They are represented as dictionaries:

.. code-block:: json

    {
        "position": int,
        "offset": int,
        "region": str
    }

Where:

- **position**: a positive integer
- **offset**: an integer indicating the offset relative to the position
- **region**: a string describing the region type (``""`` for standard coding positions, ``"-"`` for 5' UTR, ``"*"`` for 3' UTR, ``"u"`` for upstream and ``"d"`` for downstream)

Coding Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from mutalyzer_crossmapper import Coding
    >>> crossmap = Coding(_exons, _cds)
    >>> crossmap.coordinate_to_coding(31)
    {"position": -1, "offset": 0, "region": "-"}
    >>> crossmap.coding_to_coordinate({"position": -1, "offset": 0, "region": "-"})
    31

Notes
~~~~~

- The flag ``inverted=True`` can be used for transcripts on the reverse complement strand.

Protein
-------

Protein Position Model
~~~~~~~~~~~~~~~~~~~~~~

Protein positions follow the HGVS ``p`` coordinate system. They are represented as dictionaries:

.. code-block:: json

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
- **region**: a string describing the region type (``""`` for standard positions)

Protein Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Conversions between protein positions and coordinates:

.. code-block:: python

    >>> crossmap.coordinate_to_protein(41)
    {"position": 2, "position_in_codon": 2, "offset": 1, "region": ""}
    >>> crossmap.protein_to_coordinate({"position": 2, "position_in_codon": 2, "offset": 1, "region": ""})
    41


.. _numbering: http://varnomen.hgvs.org/bg-material/numbering/
.. _ReadTheDocs: https://mutalyzer-crossmapper.readthedocs.io
