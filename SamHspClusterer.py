# =========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
# Author: Siyan Liu
# =========================================================================
from __future__ import (absolute_import, division, print_function,
                        unicode_literals, generators, nested_scopes, with_statement)


class SamHspClusterer:
    """
    This class eliminates overlaps in HSPs by choosing a high-scoring subset
    with no overlaps.

    Attributes:
        none
    Instance Methods:
        clusterer=SamHspClusterer()
    Cla Methods:
        clustered=SamHspClusterer.cluster(HSPs)
    """

    def __init__(self):
        pass

    @classmethod
    def overlap(cls, hsp, HSPs):
        """
        This method determines whether a given HSP overlaps any of the
        other HSPs in a given list.
        """
        for other in HSPs:
            if hsp.overlapsOnRead(other):
                return True
        return False

    @classmethod
    def cluster(cls, raw):
        """
        This method picks a nonoverlapping set of the highest-scoring HSPs.
        """
        nonoverlapping = []
        HSPs = [x for x in raw]
        for hsp in HSPs:
            hsp.computeScore()
        HSPs.sort(key=lambda hsp: -hsp.getScore())
        for hsp in HSPs:
            if not cls.overlap(hsp, nonoverlapping):
                nonoverlapping.append(hsp)
        nonoverlapping.sort(key=lambda hsp: hsp.getReadInterval().getBegin())
        return nonoverlapping
