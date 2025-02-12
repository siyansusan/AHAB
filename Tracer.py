# =========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
# Author: Siyan Liu
# =========================================================================
from __future__ import (absolute_import, division, print_function,
                        unicode_literals, generators, nested_scopes, with_statement)

from Rex import Rex
rex = Rex()
from ConfigFile import ConfigFile
from Strand import Strand


class Tracer:
    """
    This class should be used to encapsulate common functionality across all
    different analysis scripts.

    Attributes:
        readsBinned : int
    Instance Methods:
        tracer=Tracer(OUTPUT_DIR)
        tracer.bin(Annotation,FILE)
        tracer.dump(Annotation,FILE)
        tracer.getAlignabilities(anno)
    Class Methods:
        none
    Private methods:
    """

    def __init__(self, configFile):
        self.config = ConfigFile(configFile)
        self.readsBinned = 0
        self.CHROMS = set()

    def dump(self, anno, FILE):
        """
        This method prints out debugging information for the HSPs of a read.
        """
        HSPs = anno.getHSPs()
        numHSPs = len(HSPs)
        print(anno.getReadID(), numHSPs, sep="\t", file=FILE, flush=True)
        for hsp in HSPs:
            print("\t",
                  hsp.getRefName(),
                  Strand.toString(hsp.getStrand()),
                  hsp.getReadInterval().toString(),
                  hsp.getRefInterval().toString(),
                  hsp.getCigar().toString(),
                  hsp.getPercentIdentity(),
                  hsp.getAlignability(),
                  hsp.getSeq(),
                  sep="\t", file=FILE, flush=True)

    def bin(self, anno, FILE):
        """
        This method bins a read by writing into a bin file.
        """
        readSeq = anno.getSamRecord().getSequence()
        print(anno.getReadID(), readSeq, sep="\t", file=FILE, flush=True)
        self.readsBinned += 1

    def getMinAlignability(self, A):
        """
        Find the minimum alignability.
        """
        a = []
        for x in A:
            if x is not None:
                a.append(x)
        if len(a) > 0:
            return min(a)
        return 0


    def getAlignabilities(self, anno):
        """
        This calls bigwit.stats() to get the alignabilities for all windows
        overlapping any HSPs in this annotation
        """
        for hsp in anno.getHSPs():
            if hsp.getRefName() in self.CHROMS:
                refCoords = hsp.getRefInterval()
                stats = self.bigwig.stats(hsp.getRefName(), refCoords.getBegin(),
                                          refCoords.getEnd(), type="mean")
                minValue = self.getMinAlignability(stats)
                hsp.setAlignability(minValue)
