Library
=======

The library provides a number of classes to perform the various conversions.


The ``Crossmap`` class
----------------------

.. code:: python

    >>> from crossmapper import Crossmap

.. list-table:: Constructor parameters.
   :header-rows: 1

   * - name
     - description
   * - ``locations``
     - List of locations.
   * - ``cds``
     - ...
   * - ``inverted``
     - Orientation.


No transcripts
^^^^^^^^^^^^^^

.. code:: python

    >>> crossmap = Crossmap()

``coordinate_to_genomic()`` and ``genomic_to_coordinate`` functions.

.. code:: python

    >>> crossmap.coordinate_to_genomic(0)
    1
    >>> crossmap.genomic_to_coordinate(1)
    0

Noncoding transcripts
^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    >>> exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
    >>> crossmap = Crossmap(exons)

``coordinate_to_noncoding()`` and ``noncoding_to_coordinate()`` functions.

.. list-table:: Noncoding positions.
   :header-rows: 1

   * - index
     - description
   * - 0
     - 
   * - 1
     - 

.. code:: python

    >>> crossmap.coordinate_to_noncoding(35)
    (14, 1)
    >>> crossmap.noncoding_to_coordinate((14, 1))
    35

``inverted`` parameter.

.. code:: python

    >>> crossmap = Crossmap(exons, inverted=True)
    >>> crossmap.coordinate_to_noncoding(35)
    (9, -1)
    >>> crossmap.noncoding_to_coordinate((9, -1))
    35

Coding transcripts
^^^^^^^^^^^^^^^^^^

.. code:: python

    >>> cds = (32, 43)
    >>> crossmap = Crossmap(exons, cds)

``coordinate_to_coding()`` and ``coding_to_coordinate()`` functions.

.. list-table:: Coding positions.
   :header-rows: 1

   * - index
     - description
   * - 0
     - 
   * - 1
     - 
   * - 2
     - 

.. code:: python

    >>> crossmap.coordinate_to_coding(31)
    (-1, 0, 0)
    >>> crossmap.coding_to_coordinate((-1, 0, 0))
    31

``inverted`` parameter.

.. code:: python

    >>> crossmap = Crossmap(exons, cds, inverted=True)
    >>> crossmap.coordinate_to_coding(31)
    (1, 0, 2)
    >>> crossmap.coding_to_coordinate((1, 0, 2))
    31


Basic classes
=============

The ``Locus`` class
-------------------

The ``MultiLocus`` class
-------------------
