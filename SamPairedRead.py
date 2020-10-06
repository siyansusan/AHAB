#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)

#=========================================================================
# class SamPairedRead
#
# This class represents a pair of SamRecord objects representing the
# two reads of a pair.
#
# Attributes:
#   read1, read2 : SamRecord
# Instance Methods:
#   pair=SamPairedRead(read1,read2)
#   L=pair.totalAlignedLength()
#   numMismatches=pair.numMismatches()
#   score=pair.matchProportion()
# Private Methods:
#   none
# Class Methods:
#   none
#=========================================================================
class SamPairedRead:
    """SamPairedRead"""
    def __init__(self,read1,read2):
        self.read1=read1
        self.read2=read2

    def getID(self):
        return self.read1.getID()

    # This returns the sum of the aligned lengths of the two reads; note
    # that this does NOT account for any overlap!
    def totalAlignedLength(self):
        L1=self.read1.getCigar().totalAlignmentLength()
        L2=self.read2.getCigar().totalAlignmentLength()
        return L1+L2

    # This returns the sum of indel bases in the two reads
    def countIndelBases(self):
        cigar1=self.read1.getCigar()
        cigar2=self.read2.getCigar()
        N1=cigar1.countIndelBases()
        N2=cigar2.countIndelBases()
        return N1+N2

    # This returns the sum of the count of mismatches between the two reads
    def numMismatches(self):
        mis1=self.read1.countMismatches()
        mis2=self.read2.countMismatches()
        return mis1+mis2

    # This computes a match proportion, assuming the two reads do NOT overlap!
    def matchProportion(self):
        L=self.totalAlignedLength()
        mis=self.numMismatches()
        score=float(L-mis)/float(L)
        return score

    # This scoring function is no longer used.
    def computeScore(self):
        mismatches=self.numMismatches()
        matches=self.totalAlignedLength()-mismatches
        numerator=matches
        indelBases=self.countIndelBases()
        denominator=1+mismatches+indelBases
        score=float(numerator)/float(denominator)
        return score



