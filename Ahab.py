#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
import os
from Rex import Rex
rex=Rex()
import pyBigWig
from ConfigFile import ConfigFile
from Strand import Strand

#=========================================================================
# Attributes:
#   readsBinned : int
# Instance Methods:
#   ahab=Ahab(OUTPUT_DIR)
#   ahab.bin(Annotation,FILE)
#   ahab.dump(Annotation,FILE)
#   ahab.getAlignabilities(anno)
# Class Methods:
#   none
# Private methods:
#=========================================================================
class Ahab:
    """Ahab"""
    def __init__(self,configFile):
        self.config=ConfigFile(configFile)
        alignabilityMapFile=self.config.lookupOrDie("ALIGNABILITY")
        self.bigwig=pyBigWig.open(alignabilityMapFile)
        self.readsBinned=0
        self.CHROMS=set()

    def dump(self,anno,FILE):
        HSPs=anno.getHSPs()
        numHSPs=len(HSPs)
        print(anno.getReadID(),numHSPs,sep="\t",file=FILE,flush=True)
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
                  sep="\t",file=FILE,flush=True)

    def bin(self,anno,FILE):
        readSeq=anno.getSamRecord().getSequence() ### temporary
        print(anno.getReadID(),readSeq,sep="\t",file=FILE,flush=True) ### temp

        #print(anno.getReadID(),file=FILE,flush=True)
        self.readsBinned+=1

    def getAlignabilities(self,anno):
        for hsp in anno.getHSPs():
            if(hsp.getRefName() in self.CHROMS):
                refCoords=hsp.getRefInterval()
                stats=self.bigwig.stats(hsp.getRefName(),refCoords.getBegin(),
                                        refCoords.getEnd(),type="min")
                minValue=min(stats)
                hsp.setAlignability(minValue)

                


