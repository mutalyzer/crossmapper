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

.. csv-table:: Coordinate to Genomic Position (0-4)
   :header: "Coordinate", "Position"

   0, 1
   1, 2
   2, 3
   3, 4
   4, 5
   ...

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

Here is the mapping of coordinates to noncoding positions:

.. csv-table:: Coordinate Mapping
   :header: "Coordinate", "Position", "Offset", "Region"

   0, 5, 0, u
   1, 4, 0, u
   2, 3, 0, u
   3, 2, 0, u
   4, 1, 0, u
   5, 1, 0, ""
   6, 2, 0, ""
   7, 3, 0, ""
   8, 3, 1, ""
   9, 3, 2, ""
   10, 3, 3, ""
   11, 4, -3, ""
   12, 4, -2, ""
   13, 4, -1, ""
   14, 4, 0, ""
   15, 5, 0, ""
   16, 6, 0, ""
   17, 7, 0, ""
   18, 8, 0, ""
   19, 9, 0, ""
   20, 9, 1, ""
   21, 9, 2, ""
   22, 9, 3, ""
   23, 9, 4, ""
   24, 9, 5, ""
   25, 10, -5, ""
   26, 10, -4, ""
   27, 10, -3, ""
   28, 10, -2, ""
   29, 10, -1, ""
   30, 10, 0, ""
   31, 11, 0, ""
   32, 12, 0, ""
   33, 13, 0, ""
   34, 14, 0, ""
   35, 14, 1, ""
   36, 14, 2, ""
   37, 14, 3, ""
   38, 15, -2, ""
   39, 15, -1, ""
   40, 15, 0, ""
   41, 16, 0, ""
   42, 17, 0, ""
   43, 18, 0, ""
   44, 18, 1, ""
   45, 18, 2, ""
   46, 18, 3, ""
   47, 19, -3, ""
   48, 19, -2, ""
   49, 19, -1, ""
   50, 19, 0, ""
   51, 20, 0, ""
   52, 20, 1, ""
   53, 20, 2, ""
   54, 20, 3, ""
   55, 20, 4, ""
   56, 20, 5, ""
   57, 20, 6, ""
   58, 20, 7, ""
   59, 20, 8, ""
   60, 20, 9, ""
   61, 21, -9, ""
   62, 21, -8, ""
   63, 21, -7, ""
   64, 21, -6, ""
   65, 21, -5, ""
   66, 21, -4, ""
   67, 21, -3, ""
   68, 21, -2, ""
   69, 21, -1, ""
   70, 21, 0, ""
   71, 22, 0, ""
   72, 1, 0, d
   73, 2, 0, d
   74, 3, 0, d
   75, 4, 0, d
   76, 5, 0, d
   77, 6, 0, d
   78, 7, 0, d
   79, 8, 0, d


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

Here is the mapping of coordinates to coding positions:

.. csv-table:: Coordinate Mapping
   :header: "Coordinate", "Position", "Offset", "Region"

   0, 5, 0, u
   1, 4, 0, u
   2, 3, 0, u
   3, 2, 0, u
   4, 1, 0, u
   5, 11, 0, "-"
   6, 10, 0, "-"
   7, 9, 0, "-"
   8, 9, 1, "-"
   9, 9, 2, "-"
   10, 9, 3, "-"
   11, 8, -3, "-"
   12, 8, -2, "-"
   13, 8, -1, "-"
   14, 8, 0, "-"
   15, 7, 0, "-"
   16, 6, 0, "-"
   17, 5, 0, "-"
   18, 4, 0, "-"
   19, 3, 0, "-"
   20, 3, 1, "-"
   21, 3, 2, "-"
   22, 3, 3, "-"
   23, 3, 4, "-"
   24, 3, 5, "-"
   25, 2, -5, "-"
   26, 2, -4, "-"
   27, 2, -3, "-"
   28, 2, -2, "-"
   29, 2, -1, "-"
   30, 2, 0, "-"
   31, 1, 0, "-"
   32, 1, 0, ""
   33, 2, 0, ""
   34, 3, 0, ""
   35, 3, 1, ""
   36, 3, 2, ""
   37, 3, 3, ""
   38, 4, -2, ""
   39, 4, -1, ""
   40, 4, 0, ""
   41, 5, 0, ""
   42, 6, 0, ""
   43, 1, 0, "*"
   44, 1, 1, "*"
   45, 1, 2, "*"
   46, 1, 3, "*"
   47, 2, -3, "*"
   48, 2, -2, "*"
   49, 2, -1, "*"
   50, 2, 0, "*"
   51, 3, 0, "*"
   52, 3, 1, "*"
   53, 3, 2, "*"
   54, 3, 3, "*"
   55, 3, 4, "*"
   56, 3, 5, "*"
   57, 3, 6, "*"
   58, 3, 7, "*"
   59, 3, 8, "*"
   60, 3, 9, "*"
   61, 4, -9, "*"
   62, 4, -8, "*"
   63, 4, -7, "*"
   64, 4, -6, "*"
   65, 4, -5, "*"
   66, 4, -4, "*"
   67, 4, -3, "*"
   68, 4, -2, "*"
   69, 4, -1, "*"
   70, 4, 0, "*"
   71, 5, 0, "*"
   72, 1, 0, d
   73, 2, 0, d
   74, 3, 0, d
   75, 4, 0, d
   76, 5, 0, d
   77, 6, 0, d
   78, 7, 0, d
   79, 8, 0, d

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
