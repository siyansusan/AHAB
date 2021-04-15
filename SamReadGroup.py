# =========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
# Author: Siyan Liu
# =========================================================================
from __future__ import (absolute_import, division, print_function,
                        unicode_literals, generators, nested_scopes, with_statement)


class SamReadGroup:
    """
    This class represents a group of paired reads having the same read ID.

    Attributes:
        ID : string
        readPairs : array of SamPairedRead
    Instance Methods:
        group=SamPairedReadGroup()
        group.addPair(pairedRead)
        id=group.getID()
        reads=group.getReads() : returns array of SamPairedRead
        n=group.numReadPairs()
        reads=group.getReadEnds(end1or2)
    Class Methods:
        None
    """

    def __init__(self):
        self.ID = None
        self.reads = []

    def __len__(self):
        return len(self.reads)

    def getReadEnds(self, end1or2):
        """
        This returns either all of the first reads (1) or the second reads (2)
        of the pairs, depending on the parameter "end1or2".
        """
        reads = []
        for read in self.reads:
            if end1or2 == 1:
                if read.flag_firstOfPair():
                    reads.append(read)
            else:
                if read.flag_secondOfPair():
                    reads.append(read)
        return reads

    def numRead(self):
        return len(self.reads)

    def getID(self):
        return self.ID

    def getReads(self):
        """
        This returns the read pairs.
        """
        return self.reads
