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
from SamRecord import SamRecord
from CigarString import CigarString

#=========================================================================
# Attributes:
#   keepOps : set of string
# Instance Methods:
#   factory=SamHspFactory()
#   HSPs=factory.makeHSPs(reads)
# Private Methods:
#   cigar=self.processCigar(cigar)
# Class Methods:
#   none
#=========================================================================
class SamHspFactory:
    """SamHspFactory"""
    def __init__(self):
        self.keepOps=set(["M","I","D","=","X"])

    def makeHSPs(self,reads):
        HSPs=[]
        for read in reads:
            cigar=read.getCigar()
            cigar.computeIntervals(read.getRefPos())
            cigar=self.processCigar(cigar)
            hsp=SamHSP(read,cigar)
            HSPs.append(hsp)
        return HSPs

    def processCigar(self,cigar):
        keepOps=self.keepOps
        keep=[]
        n=cigar.length()
        for i in range(n):
            op=cigar[i]
            #if(op.getOp()=="S" and i>0 and i<n-1):
            #    raise Exception("Nested soft mask detected in cigar string")
            if(op.getOp() in keepOps):
                keep.append(op)
        new=CigarString("")
        new.setOps(keep)
        return new
