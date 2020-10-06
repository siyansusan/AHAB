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
# This class eliminates overlaps in HSPs by choosing a high-scoring subset
# with no overlaps.
#
# Attributes:
#   none
# Instance Methods:
#   clusterer=SamHspClusterer()
# Cla Methods:
#   clustered=SamHspClusterer.cluster(HSPs)
#=========================================================================
class SamHspClusterer:
    """SamHspClusterer"""
    def __init__(self):
        pass

    # This method determines whether a given HSP overlaps any of the
    # other HSPs in a given list
    @classmethod
    def overlap(cls,hsp,HSPs):
        for other in HSPs:
            if(hsp.overlapsOnRead(other)): return True
        return False

    # This method picks a nonoverlapping set of the highest-scoring HSPs
    @classmethod
    def cluster(cls,raw):
        nonoverlapping=[]
        HSPs=[x for x in raw]
        for hsp in HSPs: hsp.computeScore()
        HSPs.sort(key=lambda hsp: -hsp.getScore()) # sort in reverse order
        for hsp in HSPs:
            if(not cls.overlap(hsp,nonoverlapping)):
                nonoverlapping.append(hsp)
        nonoverlapping.sort(key=lambda hsp: hsp.getReadInterval().getBegin())
        return nonoverlapping

