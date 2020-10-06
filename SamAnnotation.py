#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
from SamHSP import SamHSP
from Interval import Interval

#=========================================================================
# This class represents a set of HSPs (local alignments) for a single read.
#
# Attributes:
#   HSPs : array of SamHSP
# Instance Methods:
#   anno=SamAnnotation(HSPs) # makes shallow copy of array elements
#   ID=anno.getReadID()
#   n=anno.numHSPs()
#   HSPs=anno.getHSPs()
#   boolean=anno.allRefsSame() # are all references the same?
#   boolean=anno.allSameStrand() # do all HSPs map to same strand of refs?
#   refName=anno.firstRef() # returns reference name of first HSP
#   refNames=anno.getRefNames() # returns all reference names as a set
#   n=anno.numDifferentRefs() # returns number of different refs among HSPs
#   boolean=anno.anyRefsOverlap()
#   array=anno.getReadGaps(includeMargins=False)
#   array=anno.getReadGapLengths(includeMargins=False)
#   array=anno.getRefGaps()
#   array=anno.getRefGapLengths()
#   L=anno.getReadLength() # returns length of entire read (not just HSP)
#   N=anno.unalignedLength()
#   P=anno.unalignedProportion()
#   identity=anno.lowestPercentIdentity()
#   x=anno.getLowestAlignability()
#   SamRecord anno.getSamRecord()
# Class Methods:
#   none
#=========================================================================
class SamAnnotation:
    """SamAnnotation"""
    def __init__(self,HSPs):
        self.HSPs=[]
        for hsp in HSPs:
            self.HSPs.append(hsp)

    def __del__(self):
        for hsp in self.HSPs:
            del hsp

    # This returns the SAM record for the first HSP (note that different HSPs
    # will have different SAM records, but those records will share some info,
    # such as the read ID and read sequence)
    def getSamRecord(self):
        HSPs=self.HSPs
        if(len(HSPs)==0): raise Exception("No HSPs in Annotation")
        return HSPs[0].getRec()

    # Returns the aligned proportion, accounting for all HSPs
    def alignedProportion(self):
        L=self.getReadLength()
        N=self.alignedLength()
        return float(N)/float(L)

    # Returns the total aligned length, summed over all HSPs
    def alignedLength(self):
        # Prerequisite: HSPs are non-overlapping on the read
        aligned=0
        HSPs=self.HSPs
        for hsp in HSPs:
            aligned+=hsp.getReadInterval().getLength()
        return aligned

    def getReadID(self):
        HSPs=self.HSPs
        if(len(HSPs)==0): 
            raise Exception("No HSPs in SamAnnotation::getReadID()")
        return HSPs[0].getReadID()

    # Returns the lowest %identity across all the HSPs, for filtering purposes
    def lowestPercentIdentity(self):
        return min([x.getPercentIdentity() for x in self.HSPs])

    # Returns the lowest alignability across all the HSPs, for filtering
    def getLowestAlignability(self):
        return min([x.getAlignability() for x in self.HSPs])

    # Tells whether all HSPs had the same strand
    def allSameStrand(self):
        HSPs=self.HSPs
        n=len(HSPs)
        if(n==0): return True
        strand=HSPs[0].getStrand()
        for i in range(1,n):
            if(HSPs[i].getStrand()!=strand): return False
        return True

    # Returns the length of the full read
    def getReadLength(self):
        HSPs=self.HSPs
        n=len(HSPs)
        if(n==0): raise Exception("Don't know read length: no HSPs")
        return HSPs[0].getRec().seqLength()

    # Returns a vector of gaps *between* (not within!) HSPs
    def getReadGaps(self,includeMargins=False):
        L=self.getReadLength()
        HSPs=self.HSPs
        numHSPs=len(HSPs)
        if(numHSPs==0): return Interval(0,L)
        intervals=[]
        if(includeMargins):
            b=HSPs[0].getReadInterval().getBegin()
            if(b>0): intervals.append(Interval(0,b))
        for i in range(numHSPs-1):
            b=HSPs[i].getReadInterval().getEnd()
            e=HSPs[i+1].getReadInterval().getBegin()
            if(b<e): intervals.append(Interval(b,e))
        if(includeMargins):
            e=HSPs[numHSPs-1].getReadInterval().getEnd()
            if(e<L):
                intervals.append(Interval(e,L))
        return intervals

    # This is similar to getReadGaps(), except that the coordinates are
    # on the reference instead of the read
    def getRefGaps(self):
        HSPs=self.HSPs
        numHSPs=len(HSPs)
        if(numHSPs==0): return []
        gaps=[]
        for i in range(numHSPs-1):
            b=HSPs[i].getRefInterval().getEnd()
            e=HSPs[i+1].getRefInterval().getBegin()
            if(b<e): gaps.append(Interval(b,e))
        return gaps

    # Do any of the reference segments overlap?
    def anyRefsOverlap(self):
        HSPs=self.HSPs
        numHSPs=len(HSPs)
        for i in range(numHSPs-1):
            for j in range(i+1,numHSPs):
                if(HSPs[i].overlapsOnRef(HSPs[j])):
                    return True
        return False

    # Returns lengths of gaps between HSPs, measured on the read (not the
    # reference)
    def getReadGapLengths(self,includeMargins=False):
        intervals=self.getReadGaps(includeMargins)
        lengths=[x.getLength() for x in intervals]
        return lengths

    # Similar to getReadGapLengths(), but measured on the reference
    def getRefGapLengths(self):
        intervals=self.getRefGaps()
        lengths=[x.getLength() for x in intervals]
        return lengths
            
    # Returns the cardinality of the set of references this read is mapped to
    def numDifferentRefs(self):
        names=self.getRefNames()
        return len(names)

    # Returns the names of all reference sequences this read is mapped to in
    # this annotation
    def getRefNames(self):
        names=set()
        for hsp in self.HSPs:
            set.add(hsp.getRefName())
        return names

    # Returns the reference of the first HSP
    def firstRef(self):
        HSPs=self.HSPs
        n=len(HSPs)
        if(n==0): raise Exception("HSPs is empty in SamAnnotation")
        return HSPs[0].getRefName()

    def numHSPs(self):
        return len(self.HSPs)

    def getHSPs(self):
        return self.HSPs

    # Tells whether all HSPs in this annotation map to the same reference
    def allRefsSame(self):
        HSPs=self.HSPs
        n=len(HSPs)
        if(n==0): return True
        ref=HSPs[0].getRefName()
        for i in range(1,n):
            if(HSPs[i].getRefName()!=ref):
                return False
        return True


