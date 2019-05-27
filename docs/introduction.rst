Introduction
============

Converting between the transcript oriented c. or n. and the genomic oriented g.
numbering systems can be difficult, especially when the transcript in question
resides on the reverse complement strand.

There is no zero (0) in any of the HGVS numbering_ systems. While this is not a
major problem in itself, when combined with signed integers in the c. and n.
numbering systems (e.g., position -1 and 1 are adjacent), it gives rise to
various off by one errors when the conversion is not done properly. The use of
*offsets* for intronic positions (e.g., position 12 and 12+1 can be adjacent)
introduce yet an other type of discontinuity which makes doing arithmetical
operations in these numbering systems extremely tedious and error prone.
Finally, for transcripts that reside on the reverse complement strand, the
direction of the numbering system is opposite to that of the genomic one (e.g.,
if c.1 equals g.10, then c.2 equals g.9), which introduces yet an other level
of complexity.

This library aims to solve all aforementioned problems by providing an
interface that is able to convert from any HGVS numbering system to a
conventional *coordinate* system and back.


.. _numbering: http://varnomen.hgvs.org/bg-material/numbering/
