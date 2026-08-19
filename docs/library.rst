Library
=======

The package provides conversion helpers between zero-based genomic
coordinates and HGVS-like point models for genomic, non-coding, coding,
and protein contexts.

Coordinate And Location Conventions
-----------------------------------

- Coordinates are zero-based integers.
- Locations are provided as half-open intervals: ``(start, end)`` with
  ``start`` inclusive and ``end`` exclusive.
- HGVS-style positions exposed by public dataclasses are one-based.


The ``Genomic`` class
---------------------

The ``Genomic`` class provides an interface to conversions between genomic
(``g.``, ``m.``, ``o.``) positions and coordinates.

The ``GenomicPoint`` dataclass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Genomic positions follow the HGVS genomic coordinate system.
They are represented as 1-attribute dataclasses. Below is an example of ``g.1`` in HGVS.

.. code-block:: python

    >>> from mutalyzer_crossmapper import GenomicPoint
    >>> GenomicPoint(position=1)

Where:

- **position**: an integer representing a nucleotide position (> 0)

Genomic Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from mutalyzer_crossmapper import Genomic
    >>> crossmap = Genomic()

The functions ``coordinate_to_genomic()`` and ``genomic_to_coordinate()`` can be
used to convert to and from genomic positions.

.. code:: python

    >>> crossmap.coordinate_to_genomic(0)
    GenomicPoint(position=1)
    >>> crossmap.genomic_to_coordinate(GenomicPoint(position=1))
    0

See section :doc:`api/crossmap` for a detailed description.

The ``NonCoding`` class
-----------------------

On top of the functionality provided by the ``Genomic`` class, the
``NonCoding`` class provides an interface to conversions between noncoding
(``n.``, ``r.``) positions and coordinates. Conversions between positioning
systems should be done via a coordinate.

The ``NonCodingPoint`` Dataclass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``NonCodingPoint`` fields:

- ``position``: positive integer (1-based transcript position)
- ``offset``: integer intronic/outside offset
- ``region``: one of ``''``, ``'u'``, ``'d'``

.. code-block:: python

    >>> from mutalyzer_crossmapper import NonCodingPoint
    >>> NonCodingPoint(position=14, offset=1, region='')
    NonCodingPoint(position=14, offset=1, region='')

Non-Coding Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from mutalyzer_crossmapper import NonCoding
    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> crossmap = NonCoding(exons)
    >>> crossmap.coordinate_to_noncoding(35)
    NonCodingPoint(position=14, offset=1, region='')
    >>> crossmap.noncoding_to_coordinate(NonCodingPoint(position=14, offset=1, region=''))
    35

Upstream and downstream positions are represented using ``region='u'`` and
``region='d'``:

.. code-block:: python

    >>> crossmap.coordinate_to_noncoding(2)
    NonCodingPoint(position=1, offset=-3, region='u')
    >>> crossmap.coordinate_to_noncoding(73)
    NonCodingPoint(position=22, offset=2, region='d')

For reverse-complement orientation, set ``inverted=True``:

.. code-block:: python

    >>> reverse = NonCoding(exons, inverted=True)
    >>> reverse.coordinate_to_noncoding(35)
    NonCodingPoint(position=9, offset=-1, region='')

See :doc:`api/crossmap` for full API details.


The ``Coding`` Class
--------------------

``Coding`` extends ``NonCoding`` with coding DNA position logic
(``c.``, ``r.``), using exon locations and one CDS interval.

The ``CodingPoint`` Dataclass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CodingPoint`` fields:

- ``position``: positive integer
- ``offset``: integer
- ``region``: one of ``''``, ``'u'``, ``'d'``, ``'-'``, ``'*'``

.. code-block:: python

    >>> from mutalyzer_crossmapper import CodingPoint
    >>> CodingPoint(position=1, offset=3, region='*')
    CodingPoint(position=1, offset=3, region='*')

Coding Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from mutalyzer_crossmapper import Coding
    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> cds = (32, 43)
    >>> crossmap = Coding(exons, cds)
    >>> crossmap.coordinate_to_coding(31)
    CodingPoint(position=1, offset=0, region='-')
    >>> crossmap.coding_to_coordinate(CodingPoint(position=1, offset=0, region='-'))
    31

The optional ``degenerate=True`` argument maps outside-transcript ``u``/``d``
representations to degenerate ``-``/``*`` positions:

.. code-block:: python

    >>> crossmap.coordinate_to_coding(4)
    CodingPoint(position=11, offset=-1, region='u')
    >>> crossmap.coordinate_to_coding(4, degenerate=True)
    CodingPoint(position=12, offset=0, region='-')


Protein Conversion
------------------

``Coding`` also exposes conversion to and from protein-position models.

The ``ProteinPoint`` Dataclass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``ProteinPoint`` extends ``CodingPoint`` with:

- ``position_in_codon``: one of ``1``, ``2``, ``3``

.. code-block:: python

    >>> from mutalyzer_crossmapper import ProteinPoint
    >>> ProteinPoint(position=1, position_in_codon=3, offset=0, region='')
    ProteinPoint(position=1, offset=0, region='', position_in_codon=3)

.. code-block:: python

    >>> from mutalyzer_crossmapper import Coding
    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> cds = (32, 43)
    >>> crossmap = Coding(exons, cds)
    >>> crossmap.coordinate_to_protein(41)
    ProteinPoint(position=2, offset=0, region='', position_in_codon=2)
    >>> crossmap.protein_to_coordinate(
    ...     ProteinPoint(position=2, offset=0, region='', position_in_codon=2)
    ... )
    41

Note that canonical HGVS ``p.`` positions correspond to points with
``offset == 0`` and ``region == ''``.

See :doc:`api/crossmap` for full API details.


Location Helper
---------------

``nearest_location()`` finds the closest location index for a coordinate.
When two boundaries are equally close, ``p`` controls tie-breaking
(``0``: left, ``1``: right).

.. code-block:: python

    >>> from mutalyzer_crossmapper import nearest_location
    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> nearest_location(exons, 37)
    2
    >>> nearest_location(exons, 37, p=1)
    3

See :doc:`api/location` for full API details.


Basic Classes
-------------

These lower-level classes are used by ``NonCoding`` and ``Coding``.

The ``Point`` Dataclass
~~~~~~~~~~~~~~~~~~~~~~~

``Point`` is the internal coordinate model used by ``Locus`` and
``MultiLocus``.

- ``position``: zero-based position within a locus or concatenated loci
- ``offset``: relative offset
- ``region``: one of ``''``, ``'u'``, ``'d'`` (mainly for ``MultiLocus``)

See :doc:`api/dataclass`.


The ``Locus`` Class
~~~~~~~~~~~~~~~~~~~

``Locus`` maps one genomic interval to/from ``Point``.

.. code-block:: python

    >>> from mutalyzer_crossmapper.locus import Locus, Point, Coord
    >>> locus = Locus((10, 20))
    >>> locus.to_position(9)
    Point(position=0, offset=-1)
    >>> locus.to_coordinate(Point(position=0, offset=-1))
    9

Set ``inverted=True`` for reverse-complement orientation.

See :doc:`api/locus` for full API details.


The ``MultiLocus`` Class
~~~~~~~~~~~~~~~~~~~~~~~~

``MultiLocus`` maps coordinates across multiple intervals to/from a unified
``Point`` model.

.. code-block:: python

    >>> from mutalyzer_crossmapper import MultiLocus, Point
    >>> multilocus = MultiLocus([(10, 20), (40, 50)])
    >>> multilocus.to_position(22)
    Point(position=9, offset=3, region='')
    >>> multilocus.to_coordinate(Point(position=9, offset=3, region=''))
    22
    >>> multilocus.to_position(38)
    Point(position=10, offset=-2, region='')
    >>> multilocus.to_coordinate(Point(position=10, offset=-2, region=''))
    38

See :doc:`api/multi_locus` for full API details.
