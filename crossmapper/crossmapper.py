"""
Module for conversion from genomic coordinates to coding sequence
orientated coordinates and vice versa.
The conversions are done based upon a list of splice sites, the CDS start
and stop and the orientation of a transcript.
"""
from __future__ import division, unicode_literals
from future.builtins import str


def _plus(a, b):
    """
    Return a + b unless a is smaller than zero and the result is larger than
    zero, in that case it returns a + b + 1.

    In effect the number 0 is skipped while adding.
    """
    if a <= 0 and r >= 0:
        return a + b + 1
    return a + b


def _minus(a, b):
    """
    Return a - b unless a is larger than zero and the result is smaller than
    zero, in that case it returns (a - b) - 1.

    In effect the number 0 is skipped while subtracting.
    """
    if a >= 0 and r <= 0:
        return a - b - 1
    return a - b


def _minusr(a, b):
    """
    Return a - b unless a is larger than zero and b is smaller than zero, in
    that case it returns (a - b) - 1.

    In effect the number 0 is skipped while subtracting.
    """
    if a > 0 and b < 0:
        return a - b - 1
    return a - b


class Crossmap():
    """
    Convert from I{g.} to I{c.} or I{n.} notation or vice versa.
    """
    def __init__(self, rna, cds, orientation):
        """
        Initialise the class and do the cross mapping of the splice sites.

        :arg list rna: The list of RNA splice sites.
        :arg list cds: CDS start and stop (if present, may be empty).
        :arg int orientation: The orientation of the transcript
            (1: forward, E{-}1: reverse).
        """
        self._stop = None
        self._rna_length = len(rna)
        self._crossmapping = self._rna_length * [None]
        self.rna = list(rna)
        self.cds = list(cds)
        self.orientation = orientation

        self._crossmap_splice_sites()

        start = (orientation - 1) // 2
        self._trans_start = self._crossmapping[start]
        self._trans_end = self._crossmapping[start - self.orientation]

        if not self._stop:
            self._stop = self._trans_end

    def _crossmap_splice_sites(self):
        """
        Calculate either:
            1. The I{c.} notation of the CDS start and stop, including splice
               sites.
            2. The I{c.} notation of the RNA splice sites.
            3. The I{n.} notation of the RNA splice sites.

        For option 1, only provide an list with CDS splice sites.
        For option 2, provide an list with RNA splice sites and one with
            the CDS start and stop.
        For option 3, only provide an list with RNA splice sites.

        Examples:

        Crossmap(rna, [], 1)
            - Get the I{n.} notation of the RNA splice sites. The input is in
              forward notation.

        Crossmap(cds, [], 1)
            - Get the I{c.} notation of the CDS start and stop, and the
              internal splice sites. The input is in forward notation.

        Crossmap(rna, cds, -1)
            - Get the I{c.} notation of the RNA splice sites. The input is in
              reverse complement.


        The output is straightforward, except for the I{c.} notation of the
        downstream RNA splice sites. This is denoted by _stop + the distance to
        the stop codon, as an alternative to the *-notation.
        """
        c_pos = 1             # One unless we both have mRNA and CDS.
        d = self.orientation
        c = (d - 1) // -2     # c, x and y are used to unify forward and
        x = (-d - 1) // -2    # reverse complement.
        y = c * (self._rna_length - 1)

        if self.cds:
            # We both have mRNA and CDS, so we have to search for CDS start.
            i = y - c
            while d * (self.rna[i] - ((i + 1) % 2)) < d * self.cds[c] + c:
                i += d
            # Get the right boundary.
            c_pos = d * (self.rna[i] - self.cds[c] + (d * 2))

            # Go back to exon 1.
            while d * i > d * y:
                c_pos = _minus(c_pos, d * (self.rna[i] - self.rna[i - d] + d))
                i -= d * 2

        i = y - c
        while d * i < d * (y - c) + self._rna_length:
            self._crossmapping[i + c] = c_pos
            if i % 2:                  
                # We are skipping an intron, so only add 1 (mind the 0).
                c_pos = _plus(c_pos, 1)
            else:
                # We are skipping an exon, so add the length.
                c_pos = _plus(c_pos, self.rna[i + 1] - self.rna[i])

            # Set _stop when we find CDS stop.
            if self.cds and not self._stop and \
               d * self.rna[i - c + 1] >= d * self.cds[x]:
                self._stop = c_pos - (d * (self.rna[i - c + 1] - self.cds[x]))
            i += d

    def g2x(self, a):
        """
        This function calculates either:
            1. The I{n.} notation from a I{g.} notation.
            2. The I{c.} notation from a I{g.} notation.

        For option 1, only provide an array with mRNA splice sites and one with
            the I{c.} notation of the splice sites.
        For option 2, provide an array with mRNA splice sites, one with the
            I{c.} notation of the splice sites and an array with the CDS start
            and stop.

        Examples:

        Crossmap(rna, [], 1)
        g2x(i)
            - Get the I{n.} notation of a I{g.} position i. The input is in
              forward notation.

        Crossmap(rna, cds, -1);
        g2x(i);
            - Get the I{c.} notation of a I{g.} position i. The input is in
              reverse notation.

        The output is fully compatible with the HVGS nomenclature as defined
        on 01-07-2009.

        :arg int a: The genomic position that must be translated.

        :returns str: The I{c.} or I{n.} notation of position a.
        """
        # TODO update documentation.
        d = self.orientation
        c = (d - 1) // -2     # c and y are used to unify forward and reverse
        y = c * (self._rna_length - 1) # complement.

        if d * a < d * self.rna[y]:
            # A position before the first exon.
            return ((self._crossmapping[y]), -d * (self.rna[y] - a))
        if d * a > d * self.rna[self._rna_length - y - 1]:
            # After the last exon.
            return (self._crossmapping[self._rna_length - y - 1],
                    d * (a - self.rna[self._rna_length - y - 1]))

        for i in xrange(self._rna_length):
            # A "normal" position.
            if i % 2:
            # We're checking the intron positions.
                if self.rna[i] < a and a < self.rna[i + 1]: # Intron.
                    if d * (a - self.rna[i]) > d * (self.rna[i + 1] - a):
                        # The position was closer to the next exon.
                        return (self._crossmapping[i + 1 - c],
                                -d * (self.rna[i + 1 - c] - a))
                    # The position was closer to the previous exon.
                    return (self._crossmapping[i + c],
                            d * (a - self.rna[i + c]))
            else:
            # We're checking the exon positions.
                if self.rna[i] <= a and a <= self.rna[i + 1]:
                    return (_plus(self._crossmapping[i + c],
                                        d * (a - self.rna[i + c])), 0)

    def x2g(self, a, b):
        """
        This function calculates either:
            1. The I{g.} notation from a I{n.} notation.
            2. The I{g.} notation from a I{c.} notation.

        Whether option 1 or 2 applies depends on the content of mRNAm.

        Examples:

        Crossmap(rna, [], 1)
        x2g(i)
            - Get the I{g.} notation of a I{n.} position i. The input is in
              forward notation.

        Crossmap(rna, cds, -1);
        x2g(i, j);
            - Get the I{g.} notation of a I{c.} position i with offset j. The
              input is in reverse notation.

        :arg int a: The I{n.} or I{c.} position to be translated.
        :arg int b: The offset of position a.

        :returns int: A I{g.} position.
        """
        d = self.orientation
        c = (-d - 1) // -2 # Used to unify forward and reverse complement.

        # Assume a position before exon 1.
        ret = self.rna[0] - d * (self._crossmapping[0] - a)
        if d * a > d * self._crossmapping[self._rna_length - 1]:
            # It is after the last exon.
            ret = self.rna[self._rna_length - 1] + \
                  d * (a - self._crossmapping[self._rna_length - 1])
        for i in range(0, self._rna_length, 2):
            # Is it in an exon?
            if d * self._crossmapping[i] <= d * a and \
               d * a <= d * self._crossmapping[i + 1]:
                ret = self.rna[i + c] - d * \
                      _minusr(self._crossmapping[i + c], a)
        ret += d * b # Add the intron count.

        # Patch for CDS start on first nucleotide of exon 1.
        if a < 0 and self._crossmapping[d - c] == 1:
            ret += d

        return ret

    def int2main(self, a):
        """
        Convert the _stop notation to the '*' notation.

        :arg int a: An integer in _stop notation.

        :returns str: The converted notation (may be unaltered).
        """
        if a > self._stop:
            return '*' + str(a - self._stop)
        return str(a)

    def main2int(self, s):
        """
        Convert the '*' notation to the _stop notation.

        :arg str s: A string in '*' notation.

        :returns int: The converted notation (may be unaltered).
        """
        if s[0] == '*':
            return self._stop + int(s[1:])
        return int(s)

    def int2offset(self, t, fuzzy=False):
        """
        Convert a tuple of integers to offset-notation. This adds a `+',
        and `u' or `d' to the offset when appropriate. The main value is
        not returned.

        :arg tuple t: A tuple of integers: (main, offset) in _stop notation.
        :arg bool fuzzy: Denotes that the coordinate is fuzzy (i.e., offset is
            unknown).

        :returns str: The offset in HGVS notation.
        """
        if t[1] > 0:
            # The exon boundary is downstream.
            if fuzzy:
                return '+?'
            if t[0] >= self._trans_end:
                # It is downstream of the last exon.
                return "+d" + str(t[1])
            return '+' + str(t[1])
        if t[1] < 0:
            # The exon boundary is uptream.
            if fuzzy:
                return '-?'
            if t[0] <= self._trans_start:
                # It is upstream of the first exon.
                return "-u" + str(-t[1])
            return str(t[1])
        # No offset was given.
        return ''

    def offset2int(self, s):
        """
        Convert an offset in HGVS notation to an integer. This removes
        `+', `u' and `d' when present. It also converts a `?' to something
        sensible.

        :arg str s: An offset in HGVS notation.

        :returns int: The offset as an integer.
        """
        if not s:
            # No offset given.
            return 0
        if s == '?':
            # Here we ignore an uncertainty.
            return 0 # FIXME: this may have to be different.
        if s[1] == 'u' or s[1] == 'd':
            # Remove `u' or `d'.
            if s[0] == '-':
                # But save the `-'.
                return -int(s[2:])
            return int(s[2:])
        if s[1:] == '?':
            # Here we ignore an uncertainty in the intron.
            return 0 # FIXME: this may have to be different.
        if s[0] == '-':
            # Save the `-' here too.
            return -int(s[1:])
        return int(s[1:])

    def tuple2string(self, t, fuzzy=False):
        """
        Convert a tuple (main, offset) in _stop notation to I{c.} notation.

        :arg tuple t: A tuple (main, offset) in _stop notation.
        :arg bool fuzzy: Denotes that the coordinate is fuzzy (i.e., offset is
            unknown).

        :returns str: The position in HGVS notation.
        """
        if t[0] >= self._trans_end or t[0] <= self._trans_start:
            return str(self.int2main(_minus(t[0], -t[1])))
        return str(self.int2main(t[0])) + str(self.int2offset(t, fuzzy))

    def g2c(self, a, fuzzy=False):
        """
        Uses both g2x() and tuple2string() to translate a genomic position
        to _stop notation to I{c.} notation.

        :arg int a: The genomic position that must be translated.
        :arg bool fuzzy: Denotes that the coordinate is fuzzy (i.e., offset is
            unknown).

        :returns str: The position in HGVS notation.
        """
        return self.tuple2string(self.g2x(a), fuzzy)

    def info(self):
        """
        Return transcription start, transcription end and CDS stop.

        :returns tuple: (trans_start, trans_stop, cds_stop).
        """
        return (self._trans_start, self._trans_end, self._stop)

    def getSpliceSite(self, number):
        """
        Return the coordinate of a splice site.

        :arg int number: the number of the RNA splice site counting from
            transcription start.
        :returns int: coordinate of the RNA splice site.
        """
        if self.orientation == 1:
            return int(self.rna[number])
        return int(self.rna[len(self.rna) - number - 1])

    def numberOfIntrons(self):
        """
        Returns the number of introns.

        :returns int: number of introns.
        """
        return len(self.rna) // 2 - 1

    def numberOfExons(self):
        """
        Returns the number of exons.

        :returns int: number of exons.
        """
        return len(self.rna) // 2
