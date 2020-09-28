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
# Attributes:
#   ID : string
#   readPairs : array of SamPairedRead
# Instance Methods:
#   group=SamPairedReadGroup()
#   group.addPair(pairedRead)
#   id=group.getID()
#   reads=group.getReads() : returns array of SamPairedRead
#   n=group.numReadPairs()
#   reads=group.getReadEnds(end1or2)
# Class Methods:
#   
#=========================================================================
class SamPairedReadGroup:
    """SamPairedReadGroup"""
    def __init__(self):
        self.ID=None
        self.readPairs=[]

    def getReadEnds(self,end1or2):
        reads=[]
        for pair in self.readPairs:
            if(end1or2==1): reads.append(pair.read1)
            else: reads.append(pair.read2)
        return reads

    def numReadPairs(self):
        return len(self.readPairs)

    def getID(self):
        return self.ID

    def getReads(self):
        return self.readPairs

    def addPair(self,pair):
        if(len(self.readPairs)==0):
            self.ID=pair.read1.getID()
        self.readPairs.append(pair)


