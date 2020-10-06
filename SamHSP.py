#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
from CigarString import CigarString
from Interval import Interval
from Strand import Strand

#=========================================================================
# Class SamHSP
#
# This class represents an HSP (High-scoring Segment Pair), which is a local
# alignment between a read and a reference sequence.  It encapsulates the
# coordinates on the read and reference, and information about the quality of
# that local alignment.  It also has a link back to the SamRecord used to
# make the HSP.
#
# Attributes:
#   cigar : CigarString
#   refName : string
#   strand : Strand
#   readInterval : Interval
#   refInterval : Interval
#   score : float = #matches/(1+#mismatches+#indelbases)
#   percentIdentity : float = #matches/(#matches+#mismatches+#indelbases)
#   alignability : float (from ENCODE alignability map)
#   rec : the SamRecord this HSP came from
# Instance Methods:
#   hsp=SamHSP(rec,cigar) # rec is a SamRecord
#   cigar=hsp.getCigar() # returns CigarString object
#   refName=hsp.getRefName()
#   Strand hsp.getStrand()
#   seq=hsp.getSeq()
#   boolean hsp.forwardStrand()
#   boolean=hsp.overlapsOnRead(otherHSP)
#   boolean=hsp.overlapsOnRef(otherHSP)
#   interval=hsp.getReadInterval()  # returns Interval object
#   interval=hsp.getRefInterval() # returns Interval object
#   hsp.computeScore()
#   score=hsp.getScore()
#   identity=hsp.getPercentIdentity()
#   hsp.setAlignability(x)
#   x=hsp.getAlignability()
#   str=hsp.toString()
#   rec=hsp.getRec() # returns the SamRecord this HSP came from
#   ID=hsp.getReadID()
# Private Methods:
#   self.computeIntervals()
# Class Methods:
#   none
#=========================================================================
class SamHSP:
    """SamHSP"""
    def __init__(self,rec,cigar):
        self.cigar=cigar
        self.refName=rec.getRefName()
        self.rec=rec
        self.computeIntervals()
        self.score=None
        self.percentIdentity=None
        self.strand=Strand.REVERSE if rec.flag_revComp() else Strand.FORWARD
        self.alignability=None

    def containsIndels(self):
        return not self.cigar.completeMatch()

    def getReadID(self):
        return self.rec.getID()

    def setAlignability(self,x):
        self.alignability=x

    def getAlignability(self):
        return self.alignability

    def getSeq(self):
        return self.rec.getSequence()[self.readInterval.getBegin():
                                          self.readInterval.getEnd()]

    # This returns the original SamRecord used to make this HSP
    def getRec(self):
        return self.rec

    def getStrand(self):
        return self.strand

    def forwardStrand(self):
        return self.strand==Strand.FORWARD

    # This generates a printable string for debugging
    def toString(self):
        return self.refName+"|"+Strand.toString(self.strand)+"|"+\
            self.refInterval.toString()+"|"+\
            self.readInterval.toString()+"|"+\
            self.cigar.toString()+"|"+\
            str(round(self.getPercentIdentity(),3))+"|"+\
            self.getSeq()

    def getPercentIdentity(self):
        if(self.percentIdentity is None):
            mismatches=self.rec.countMismatches()
            matches=self.cigar.totalAlignmentLength()-mismatches
            indelBases=self.cigar.countIndelBases()
            numerator=matches
            denominator=matches+mismatches+indelBases
            self.percentIdentity=float(numerator)/float(denominator)
        return self.percentIdentity

    def getScore(self):
        return self.score

    # This computes a score that is different from the percent identity:
    # it is #matches / (1+#mismatches+#indelbases).  This score is used
    # in the clustering step to eliminate overlapping HSPs of lower quality.
    def computeScore(self):
        mismatches=self.rec.countMismatches()
        matches=self.cigar.totalAlignmentLength()-mismatches
        indelBases=self.cigar.countIndelBases()
        numerator=matches
        denominator=1+mismatches+indelBases
        self.score=round(float(numerator)/float(denominator),2)

    # Do these HSPs overlap on the read?
    def overlapsOnRead(self,other):
        return self.readInterval.overlaps(other.readInterval)

    # Do these HSPs overlap on the reference sequence?
    def overlapsOnRef(self,other):
        return self.refInterval.overlaps(other.refInterval)

    # This returns the *processed* CIGAR string, which does not include the
    # soft-mask elements
    def getCigar(self):
        return self.cigar

    def getRefName(self):
        return self.refName

    def getReadInterval(self):
        return self.readInterval

    def getRefInterval(self):
        return self.refInterval

    def computeIntervals(self):
        cigar=self.cigar
        n=cigar.length()
        firstOp=cigar[0]
        lastOp=cigar[n-1]
        self.readInterval=Interval(firstOp.getQueryInterval().getBegin(),
                                   lastOp.getQueryInterval().getEnd())
        self.refInterval=Interval(firstOp.getRefInterval().getBegin(),
                                  lastOp.getRefInterval().getEnd())
            


