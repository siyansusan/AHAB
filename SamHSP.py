# =========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
# Author: Siyan Liu
# =========================================================================
from __future__ import (absolute_import, division, print_function,
                        unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
                      chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)

from CigarString import CigarString

from Interval import Interval
from Strand import Strand


# =========================================================================


class SamHSP:
    """
    This class represents an HSP (High-scoring Segment Pair), which is a local
    alignment between a read and a reference sequence.  It encapsulates the
    coordinates on the read and reference, and information about the quality of
    that local alignment.  It also has a link back to the SamRecord used to
    make the HSP.

    Attributes:
        cigar : CigarString
        refName : string
        strand : Strand
        readInterval : Interval
        refInterval : Interval
        score : float = #matches/(1+#mismatches+#indelbases)
        percentIdentity : float = #matches/(#matches+#mismatches+#indelbases)
        alignability : float (from ENCODE alignability map)
        rec : the SamRecord this HSP came from
    Instance Methods:
        hsp=SamHSP(rec,cigar) # rec is a SamRecord
        cigar=hsp.getCigar() # returns CigarString object
        refName=hsp.getRefName()
        Strand hsp.getStrand()
        seq=hsp.getSeq()
        boolean hsp.forwardStrand()
        boolean=hsp.overlapsOnRead(otherHSP)
        boolean=hsp.overlapsOnRef(otherHSP)
        interval=hsp.getReadInterval()  # returns Interval object
        interval=hsp.getRefInterval() # returns Interval object
        hsp.computeScore()
        score=hsp.getScore()
        identity=hsp.getPercentIdentity()
        hsp.setAlignability(x)
        x=hsp.getAlignability()
        str=hsp.toString()
        rec=hsp.getRec() # returns the SamRecord this HSP came from
        ID=hsp.getReadID()
    Private Methods:
        self.computeIntervals()
    Class Methods:
        none
    """

    def __init__(self, rec, cigar):
        self.cigar = cigar
        self.refName = rec.getRefName()
        self.rec = rec
        self.computeIntervals()
        self.score = None
        self.percentIdentity = None
        self.strand = Strand.REVERSE if rec.flag_revComp() else Strand.FORWARD
        self.alignability = None

    def containsOnTargetIndels(self, cut_site):
        """
        This method returns true if HSP contains indel 15+- from the cut site.
        """
        begin_pos = int(self.refInterval.getBegin())
        for op in self.cigar.ops:
            if op.getOp() in ("I", "D"):
                op_interval = Interval(begin_pos, begin_pos + op.getLength())
                for site in cut_site:
                    cut_range = Interval(int(site) - 15, int(site + 15))
                    if cut_range.overlaps(op_interval):
                        return True
            begin_pos = begin_pos + op.getLength()
        return False

    def containsIndels(self):
        return not self.cigar.completeMatch()

    def getReadID(self):
        return self.rec.getID()

    def setAlignability(self, x):
        self.alignability = x

    def getAlignability(self):
        return self.alignability

    def getSeq(self):
        return self.rec.getSequence()[self.readInterval.getBegin():
                                      self.readInterval.getEnd()]

    def getRec(self):
        """
        This returns the original SamRecord used to make this HSP
        """
        return self.rec

    def getStrand(self):
        return self.strand

    def forwardStrand(self):
        return self.strand == Strand.FORWARD

    def toString(self):
        """
        This generates a printable string for debugging.
        """
        return self.refName + "|" + Strand.toString(self.strand) + "|" + \
               self.refInterval.toString() + "|" + \
               self.readInterval.toString() + "|" + \
               self.cigar.toString() + "|" + \
               str(round(self.getPercentIdentity(), 3)) + "|" + \
               self.getSeq()

    def getPercentIdentity(self):
        """
        Method returns the percent of matched base pairs out of all base pairs.
        """
        if self.percentIdentity is None:
            mismatches = self.rec.countMismatches()
            matches = self.cigar.totalAlignmentLength() - mismatches
            indelBases = self.cigar.countIndelBases()
            numerator = matches
            denominator = matches + mismatches + indelBases
            self.percentIdentity = float(numerator) / float(denominator)
        return self.percentIdentity

    def getScore(self):
        return self.score

    def computeScore(self):
        """
        This computes a score that is different from the percent identity:
        it is #matches / (1+#mismatches+#indelbases).  This score is used
        in the clustering step to eliminate overlapping HSPs of lower quality.
        """
        mismatches = self.rec.countMismatches()
        matches = self.cigar.totalAlignmentLength() - mismatches
        indelBases = self.cigar.countIndelBases()
        numerator = matches
        denominator = 1 + mismatches + indelBases
        self.score = round(float(numerator) / float(denominator), 2)

    def overlapsOnRead(self, other):
        """
        Do these HSPs overlap on the read?
        """
        return (self.readInterval.begin + 10) < other.readInterval.end and (
                    other.readInterval.begin + 10) < self.readInterval.end

    def overlapsOnRef(self, other):
        """
        Do these HSPs overlap on the reference sequence?
        """
        return self.refInterval.overlaps(other.refInterval)

    def getCigar(self):
        """
        This returns the *processed* CIGAR string, which does not include the
        soft-mask elements
        """
        return self.cigar

    def getRefName(self):
        return self.refName

    def getReadInterval(self):
        return self.readInterval

    def getRefInterval(self):
        return self.refInterval

    def computeIntervals(self):
        """
        Begin and end of an HSP.
        """
        cigar = self.cigar
        n = cigar.length()
        firstOp = cigar[0]
        lastOp = cigar[n - 1]
        self.readInterval = Interval(firstOp.getQueryInterval().getBegin(),
                                     lastOp.getQueryInterval().getEnd())
        self.refInterval = Interval(firstOp.getRefInterval().getBegin(),
                                    lastOp.getRefInterval().getEnd())
