Library
=======

The library provides a number of classes to perform various conversions.


The ``Crossmap`` class
----------------------

The ``Crossmap`` class provides an interface to all conversions between
positioning systems and coordinates. Conversions between positioning systems
should be done via a coordinate.

.. code:: python

    >>> from mutalyzer_crossmapper import Crossmap

The constructor takes several optional parameters, these determine which types
of conversions can be done.

.. list-table:: Constructor parameters.
   :header-rows: 1

   * - name
     - description
     - default
   * - ``locations``
     - List of locations ``(start, end)`` of exons.
     - ``None``
   * - ``cds``
     - Location ``(start, end)`` of the CDS.
     - ``None``
   * - ``inverted``
     - Orientation.
     - ``False``

No transcripts
^^^^^^^^^^^^^^

When no parameters are provided to the constructor, only conversions concerning
genomic positions are availabe.

.. code:: python

    >>> crossmap = Crossmap()

The functions ``coordinate_to_genomic()`` and ``genomic_to_coordinate`` can be
used to convert to and from genomic positions.

.. code:: python

    >>> crossmap.coordinate_to_genomic(0)
    1
    >>> crossmap.genomic_to_coordinate(1)
    0

Noncoding transcripts
^^^^^^^^^^^^^^^^^^^^^

When a list of exon boundary coordinates is passed to the constructor,
conversions concerning noncoding transcripts are available.

.. code:: python

    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> crossmap = Crossmap(exons)

Now the functions ``coordinate_to_noncoding()`` and
``noncoding_to_coordinate()`` can be used. These functions use a 2-tuple to
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

In our example, the HGVS position "g.36" (coordinate ``35``) is equivalent to
position "n.14+1". We can convert between these two as follows.

.. code:: python

    >>> crossmap.coordinate_to_noncoding(35)
    (14, 1)
    >>> crossmap.noncoding_to_coordinate((14, 1))
    35

For transcripts that reside on the reverse complement strand, the ``inverted``
parameter should be set to ``True``. In our example, HGVS position "g.36"
(coordinate ``35``) is now equivalent to position "n.9-1".

.. code:: python

    >>> crossmap = Crossmap(exons, inverted=True)
    >>> crossmap.coordinate_to_noncoding(35)
    (9, -1)
    >>> crossmap.noncoding_to_coordinate((9, -1))
    35

Coding transcripts
^^^^^^^^^^^^^^^^^^

When both a list of exon boundary coordinates, as well as the CDS coordinates
are passed to the constructor, conversions concerning coding transcripts are
available.

.. code:: python

    >>> cds = (32, 43)
    >>> crossmap = Crossmap(exons, cds)

Now the functions ``coordinate_to_coding()`` and ``coding_to_coordinate()`` can
be used. These functions use a 3-tuple to represent a coding position.

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

The region denotes the location of the position with respect to the CDS. This
is needed in order to work with the HGVS "-" and "*" positions.

.. list-table:: Coding position regions.
   :header-rows: 1

   * - value
     - description
     - HGVS example
   * - ``0``
     - Upstream of the CDS.
     - "c.-10"
   * - ``1``
     - In the CDS.
     - "c.1"
   * - ``2``
     - Downstream of the CDS.
     - "c.*10"

In our example, the HGVS position "g.32" (coordinate ``31``) is equivalent to
position "c.-1". We can convert between these two as follows.

.. code:: python

    >>> crossmap.coordinate_to_coding(31)
    (-1, 0, 0)
    >>> crossmap.coding_to_coordinate((-1, 0, 0))
    31

Additionally, the functions ``coordinate_to_protein()`` and
``protein_to_coordinate()`` can be used. These functions use a 4-tuple to
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

In our example the HGVS position "g.42" (coordinate ``41``) corresponds with
position "p.2". We can convert between these to as follows.

.. code:: python

    >>> crossmap.coordinate_to_protein(41)
    (2, 2, 0, 1)
    >>> crossmap.protein_to_coordinate((2, 2, 0, 1))
    41

Note that the protein position only corresponds with the HGVS "p." notation
when the offset equals ``0`` and the region equals ``1``.

.. list-table:: Protein positions examples.
   :header-rows: 1

   * - coordinate
     - protein position
     - description
     - HGVS position
   * - ``31``
     - ``(-1, 3, 0, 0)``
     - Upstream position.
     - invalid
   * - ``36``
     - ``(1, 3, 2, 1)``
     - Intronic position.
     - invalid
   * - ``40``
     - ``(2, 1, 0, 1)``
     - Second amino acid, first nucleotide.
     - "p.2"
   * - ``41``
     - ``(2, 2, 0, 1)``
     - Second amino acid, second nucleotide.
     - "p.2"
   * - ``43``
     - ``(1, 1, 0, 2)``
     - Downstream position.
     - invalid


Basic classes
-------------

The ``Crossmap`` class makes use of a number of basic classes described in this
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
functions work with a 2-tuple, see the section about `Noncoding transcripts`_
for the semantics.

.. code:: python

    >>> locus.to_position(9)
    (1, -1)

For loci that reside on the reverse complement strand, the optional
``inverted`` constructor parameter should be set to ``True``.

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
