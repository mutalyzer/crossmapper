Library
=======

The library provides a number of classes to perform various conversions.


The ``Genomic`` class
---------------------

The ``Genomic`` class provides an interface to conversions between genomic
positions and coordinates.

.. code:: python

    >>> from mutalyzer_crossmapper import Genomic
    >>> crossmap = Genomic()

The functions ``coordinate_to_genomic()`` and ``genomic_to_coordinate`` can be
used to convert to and from genomic positions.

.. code:: python

    >>> crossmap.coordinate_to_genomic(0)
    1
    >>> crossmap.genomic_to_coordinate(1)
    0

See section :doc:`api/crossmap` for a detailed description.

The ``NonCoding`` class
-----------------------

On top of the functionality provided by the ``Genomic`` class, the
``NonCoding`` class provides an interface to conversions between noncoding
positions and coordinates. Conversions between positioning systems should be
done via a coordinate.

.. code:: python

    >>> from mutalyzer_crossmapper import NonCoding
    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> crossmap = NonCoding(exons)

Now the functions ``coordinate_to_noncoding()`` and
``noncoding_to_coordinate()`` can be used. These functions use a 3-tuple to
represent a noncoding position.

.. _table_noncoding:
.. list-table:: Noncoding positions.
   :header-rows: 1

   * - index
     - description
   * - 0
     - Transcript position.
   * - 1
     - Offset.
   * - 2
     - Upstream or downstream offset.

In our example, the HGVS position "g.36" (coordinate ``35``) is equivalent to
position "n.14+1". We can convert between these two as follows.

.. code:: python

    >>> crossmap.coordinate_to_noncoding(35)
    (14, 1, 0)

When the coordinate is upstream or downstream of the transcript, the last
element of the tuple denotes the offset with respect to the transcript. This
makes it possible to distinguish between intronic positions and those outside
of the transcript.

.. code:: python

    >>> crossmap.coordinate_to_noncoding(2)
    (1, -3, -3)
    >>> crossmap.coordinate_to_noncoding(73)
    (22, 2, 2)

Note that this last element is optional (and ignored) when a conversion to a
coordinate is requested.

    >>> crossmap.noncoding_to_coordinate((14, 1))
    35

For transcripts that reside on the reverse complement strand, the ``inverted``
parameter should be set to ``True``. In our example, HGVS position "g.36"
(coordinate ``35``) is now equivalent to position "n.9-1".

.. code:: python

    >>> crossmap = NonCoding(exons, inverted=True)
    >>> crossmap.coordinate_to_noncoding(35)
    (9, -1, 0)
    >>> crossmap.noncoding_to_coordinate((9, -1))
    35

See section :doc:`api/crossmap` for a detailed description.

The ``Coding`` class
--------------------

The ``Coding`` class provides an interface to all conversions between
positioning systems and coordinates. Conversions between positioning systems
should be done via a coordinate.

.. code:: python

    >>> from mutalyzer_crossmapper import Coding
    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> cds = (32, 43)
    >>> crossmap = Coding(exons, cds)

On top of the functionality provided by the ``NonCoding`` class, the functions
``coordinate_to_coding()`` and ``coding_to_coordinate()`` can be used. These
functions use a 4-tuple to represent a coding position.

.. list-table:: Coding positions.
   :header-rows: 1

   * - index
     - description
   * - 0
     - Transcript position.
   * - 1
     - Offset.
   * - 2
     - Region.
   * - 3
     - Upstream or downstream offset.

The region denotes the location of the position with respect to the CDS. This
is needed in order to work with the HGVS "-" and "*" positions.

.. list-table:: Coding position regions.
   :header-rows: 1

   * - value
     - description
     - HGVS example
   * - ``-1``
     - Upstream of the CDS.
     - "c.-10"
   * - ``0``
     - In the CDS.
     - "c.1"
   * - ``1``
     - Downstream of the CDS.
     - "c.*10"

In our example, the HGVS position "g.32" (coordinate ``31``) is equivalent to
position "c.-1". We can convert between these two as follows.

.. code:: python

    >>> crossmap.coordinate_to_coding(31)
    (-1, 0, -1, 0)
    >>> crossmap.coding_to_coordinate((-1, 0, -1))
    31

The ``coordinate_to_coding()`` function accepts an optional ``degenerate``
argument. When set to ``True``, positions outside of the transcript are no
longer described using the offset notation.

.. code:: python

    >>> crossmap.coordinate_to_coding(4)
    (-11, -1, -1, -1)
    >>> crossmap.coordinate_to_coding(4, True)
    (-12, 0, -1, -1)

Additionally, the functions ``coordinate_to_protein()`` and
``protein_to_coordinate()`` can be used. These functions use a 5-tuple to
represent a protein position.

.. list-table:: Protein positions.
   :header-rows: 1

   * - index
     - description
   * - 0
     - Protein position.
   * - 1
     - Codon position.
   * - 2
     - Offset.
   * - 3
     - Region.
   * - 4
     - Upstream or downstream offset.

In our example the HGVS position "g.42" (coordinate ``41``) corresponds with
position "p.2". We can convert between these to as follows.

.. code:: python

    >>> crossmap.coordinate_to_protein(41)
    (2, 2, 0, 0, 0)
    >>> crossmap.protein_to_coordinate((2, 2, 0, 0))
    41

Note that the protein position only corresponds with the HGVS "p." notation
when the offset equals ``0`` and the region equals ``1``. In the following
table, we show a number of annotated examples.

.. list-table:: Protein positions examples.
   :header-rows: 1

   * - coordinate
     - protein position
     - description
     - HGVS position
   * - ``4``
     - ``(-4, 2, -1, -1, -1)``
     - Upstream position.
     - invalid
   * - ``31``
     - ``(-1, 3, 0, -1, 0)``
     - 5' UTR position.
     - invalid
   * - ``36``
     - ``(1, 3, 2, 0, 0)``
     - Intronic position.
     - invalid
   * - ``40``
     - ``(2, 1, 0, 0, 0)``
     - Second amino acid, first nucleotide.
     - "p.2"
   * - ``41``
     - ``(2, 2, 0, 0, 0)``
     - Second amino acid, second nucleotide.
     - "p.2"
   * - ``43``
     - ``(1, 1, 0, 1, 0)``
     - 3' UTR position.
     - invalid
   * - ``43``
     - ``(2, 2, 2, 1, 2)``
     - Downstream position.
     - invalid

See section :doc:`api/crossmap` for a detailed description.

Locations
---------

In many cases we need to know the nearest location with respect to a
coordinate. For example, we need to know where the nearest exon is when we want
to describe a position in an intron. The ``nearest_location()`` can be used to
do exactly this.

.. code:: python

    >>> from mutalyzer_crossmapper import nearest_location
    >>> nearest_location(exons, 37)
    2
    >>> nearest_location(exons, 38)
    3

Notice that coordinate ``37`` is in the center of intron 2. By default
``nearest_location()`` will return the left location in case of a draw. This
behaviour can be altered by setting the optional argument ``p`` to ``1``.

.. code:: python

    >>> nearest_location(exons, 37, 1)
    3

See section :doc:`api/location` for a detailed description.

Basic classes
-------------

The ``Coding`` class makes use of a number of basic classes described in this
section.

The ``Locus`` class
^^^^^^^^^^^^^^^^^^^

The ``Locus`` class is used to deal with offsets with respect to a single
locus. 

.. code:: python

    >>> from mutalyzer_crossmapper import Locus
    >>> locus = Locus((10, 20))

This class provides the functions ``to_position()`` and ``to_coordinate()`` for
converting from a locus position to a coordinate and vice versa. These
functions work with a 2-tuple, see the section about `The NonCoding class`_
for the semantics.

.. code:: python

    >>> locus.to_position(9)
    (1, -1)

For loci that reside on the reverse complement strand, the optional
``inverted`` constructor parameter should be set to ``True``.

See section :doc:`api/locus` for a detailed description.

The ``MultiLocus`` class
^^^^^^^^^^^^^^^^^^^^^^^^

The ``MultiLocus`` class is used to deal with offsets with respect to multiple
loci.

.. code:: python

    >>> from mutalyzer_crossmapper import MultiLocus
    >>> multilocus = MultiLocus([(10, 20), (40, 50)])

The interface to this class is similar to that of the ``Locus`` class.

.. code:: python

    >>> multilocus.to_position(22)
    (10, 3)
    >>> multilocus.to_position(38)
    (11, -2)

See section :doc:`api/multi_locus` for a detailed description.
