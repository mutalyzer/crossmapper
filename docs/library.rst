Library
=======

The library provides a number of classes to perform various conversions.


The ``Genomic`` class
---------------------

The ``Genomic`` class provides an interface to conversions between genomic
(``g.``, ``m``, ``n``) positions and coordinates.

Genomic Position Model
~~~~~~~~~~~~~~~~~~~~~~~

Genomic positions follow the HGVS genomic coordinate system.
They are represented as 1-key dictionaries. Below is an example of `g.1` in HGVS.

.. code-block:: python

    {'position':1}

Where:

- **position**: a positive integer(>0)

Genomic Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    >>> from mutalyzer_crossmapper import Genomic
    >>> crossmap = Genomic()

The functions ``coordinate_to_genomic()`` and ``genomic_to_coordinate`` can be
used to convert to and from genomic positions.

.. code:: python

    >>> crossmap.coordinate_to_genomic(0)
    1
    >>> crossmap.genomic_to_coordinate({'position':1})
    0

See section :doc:`api/crossmap` for a detailed description.

The ``NonCoding`` class
-----------------------

On top of the functionality provided by the ``Genomic`` class, the
``NonCoding`` class provides an interface to conversions between noncoding
(``n.``, ``r.``) positions and coordinates. Conversions between positioning
systems should be done via a coordinate.

NonCoding Position Model
~~~~~~~~~~~~~~~~~~~~~~~~

Noncoding positions follow the HGVS ``n`` coordinate system. They are represented
as 3-key dictionaries. Below is an example of ``n.14+1`` in HGVS.

.. code-block:: python

    {
        'position': 14,
        'offset': 1,
        'region': ''
    }

Where:

- **position**: an interger representing a transcript position (>0)
- **offset**: an integer indicating the offset relative to the position (negative for upstream,
  positive for downstream)
- **region**: a string describing the region type (``''`` for standard, ``'u'`` for upstream,
  ``'d'`` for downstream)

NonCoding Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    >>> from mutalyzer_crossmapper import NonCoding
    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> crossmap = NonCoding(exons)

Now the functions ``coordinate_to_noncoding()`` and ``noncoding_to_coordinate()``
can be used.

In our example, the HGVS position "g.36" (coordinate ``35``) is equivalent to
position "n.14+1". We can convert between these two as follows.

.. code:: python

    >>> crossmap.coordinate_to_noncoding(35)
    {'position':14, 'offset':1, 'region':''}
    >>> crossmap.noncoding_to_coordinate({'position':14, 'offset':1, 'region':''})
    {'position':14, 'offset':1, 'region':''}

When the coordinate is upstream or downstream of the transcript,  we use ``'u`` to
present upstream and ``'d'`` to present downstream.

.. code:: python

    >>> crossmap.coordinate_to_noncoding(2)
    {'position':3, 'offset':0, 'region':'u'}
    >>> crossmap.noncoding_to_coordinate({'position':3, 'offset':0, 'region':'u'})
    2
    >>> crossmap.coordinate_to_noncoding(73)
    {'position':2, 'offset':0, 'region':'d'}
    >>> crossmap.noncoding_to_coordinate({'position':2, 'offset':0, 'region':'d'})
    73


For transcripts that reside on the reverse complement strand, the ``inverted``
parameter should be set to ``True``. In our example, HGVS position "g.36"
(coordinate ``35``) is now equivalent to position "n.9-1".

.. code:: python

    >>> crossmap = NonCoding(exons, inverted=True)
    >>> crossmap.coordinate_to_noncoding(35)
    {'position':9, 'offset':-1, 'region':''}
    >>> crossmap.noncoding_to_coordinate({'position':9, 'offset':-1, 'region':''})
    35

In the following table, we show a number of annotated examples.
.. csv-table::
   :class: table-scroll
   :header: "Coordinate", "Position", "Offset", "Region", "HGVS"

   "0", "5", "0", `u`, `n.u5`
   "4", "1", "0", `u`, `n.u1`
   "5", "1", "0", `""`, `n.1`
   "24", "9", "5", `""`, `n.9+5`
   "25", "10", "-5", `""`, `n.10-5`
   "71", "22", "0", `""`, `n.22`
   "72", "1", "0", `d`, `n.d1`
   "79", "8", "0", `d`, `n.d8`

See section :doc:`api/crossmap` for a detailed description.

The ``Coding`` class
--------------------

The ``Coding`` class provides an interface to all conversions between
coding (``c.``, ``r.``) rpositioning systems and coordinates. Conversions between
positioning systems should be done via a coordinate.

Coding Position Model
~~~~~~~~~~~~~~~~~~~~~
Coding positions follow the HGVS ``c`` coordinate system. They are
represented as 3-key dictionaries. Here is an example of ``c.*1+3``.

.. code-block:: python

    {
        'position': 1,
        'offset': 3,
        'region': '*'
    }

Where:

- **position**: an interger representing a transcript position (>0)
- **offset**: an integer indicating the offset relative to the position
- **region**: a string describing the region type (`''` for standard coding positions,
  `'-'` for 5' UTR, `'*'` for 3' UTR, `'u'` for upstream and ``'d'`` for downstream)

Coding Position Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    >>> from mutalyzer_crossmapper import Coding
    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> cds = (32, 43)
    >>> crossmap = Coding(exons, cds)

On top of the functionality provided by the ``NonCoding`` class, the functions
``coordinate_to_coding()`` and ``coding_to_coordinate()`` can be used. These
functions use a 3-key dictionary to represent a coding position.

In our example, the HGVS position "g.32" (coordinate ``31``) is equivalent to
position "c.-1". We can convert between these two as follows.

.. code:: python

    >>> crossmap.coordinate_to_coding(31)
    {'position':1, 'offset':0, 'region':'-'}
    >>> crossmap.coding_to_coordinate({'position':1, 'offset':0, 'region':'-'})
    31

The ``coordinate_to_coding()`` function accepts an optional ``degenerate``
argument. When set to ``True``, positions outside of the transcript are no
longer described using the ``'u'`` or ``'d'`` notation.

.. code:: python

    >>> crossmap.coordinate_to_coding(4)
    {'position':1, 'offset':0, 'region':'u'}
    >>> crossmap.coordinate_to_coding(4, True)
    {'position':12, 'offset':0, 'region':'-'}

In the following table, we show a number of annotated examples.

.. csv-table::
   :class:
   :header: "Coordinate", "Position", "Offset", "Region", "HGVS"

   "0", "5", "0", `u`, `c.u5`
   "4", "1", "0", `u`, `c.u1`
   "5", "11", "0", `-`, `c.-11`
   "24", "3", "5", `-`, `c.-3+5`
   "31", "1", "0", `-`, `c.-1`
   "32", "1", "0", `""`, `c.1`
   "37", "3", "3", `""`, `c.3+3`
   "38", "4", "-2", `""`, `c.4-2`
   "43", "1", "0", `*`, `c.*1`
   "61", "4", "-9", `*`, `c.*4+9`
   "71", "5", "0", `*`, `c.*5`
   "79", "8", "0", `d`, `c.d8`


Protein
-------

Additionally, the functions ``coordinate_to_protein()`` and
``protein_to_coordinate()`` can be used. These functions use a 4-key dictionary
to represent a protein position. Here is an example of ``p.1`` in HGVS.

.. code-block:: python

    {
        'position': 1,
        'position_in_codon': 3,
        'offset': 3,
        'region': ''
    }

Where:

- **position**: an interger representing the protein position (>0)
- **position_in_codon**: an integer indicating the nucleotide index within the codon (1, 2, or 3)
- **offset**: an integer indicating offset relative to the codon
- **region**: a string describing the region type (``''`` for standard positions)


In our example the HGVS position "g.42" (coordinate ``41``) corresponds with
position "p.2". We can convert between these to as follows.

.. code:: python

    >>> crossmap.coordinate_to_protein(41)
    {'position':2, 'position_in_codon':2, 'offset':0, 'region':''}
    >>> crossmap.protein_to_coordinate({'position':2, 'position_in_codon':2, offset':0, 'region':''})
    41

Note that the protein position only corresponds with the HGVS "p." notation
when the offset equals ``0`` and the region equals ``1``. In the following
table, we show a number of annotated examples.

.. csv-table::
   :class: table-scroll
   :header: "Coordinate", "Position", "position_in_codon", "Offset", "Region", "HGVS"

   "0", "4", "2", "0", `u`, ``
   "4", "4", "2", "0", `u`, ``
   "5", "4", "2", "0", `-`, ``
   "6", "4", "3", "0", `-`, ``
   "7", "3", "1", "0", `-`, ``
   "31", "1", "3", "0", `-`, ``
   "32", "1", "1", "0", ``, `p.1`
   "33", "1", "2", "0", ``, `p.1`
   "42", "2", "3", "0", ``, `p.2`
   "43", "1", "1", "0", `*`, ``
   "44", "1", "1", "1", `*`, ``
   "79", "2", "2", "0", `d`, ``



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
functions work with a 2-key dictionary, see the section about `The NonCoding class`_
for the semantics.

.. code:: python

    >>> locus.to_position(9)
    {'position':1, 'offset':-1}

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
    {'position':10, 'offset':3, 'region':''}
    >>> multilocus.to_position(38)
    {'position':11, 'offset':-2, 'region':''}

See section :doc:`api/multi_locus` for a detailed description.
