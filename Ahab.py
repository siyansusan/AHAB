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

#=========================================================================
# Attributes:
#   readsBinned : int
# Instance Methods:
#   ahab=Ahab(OUTPUT_DIR)
#   ahab.bin(readID,FILE)
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

    def bin(self,readID,FILE):
        print(readID,file=FILE,flush=True)
        self.readsBinned+=1

    def getAlignabilities(self,anno):
        for hsp in anno.getHSPs():
            refCoords=hsp.getRefInterval()
            stats=self.bigwig.stats(hsp.getRefName(),refCoords.getBegin(),
                                    refCoords.getEnd(),type="min")
            minValue=min(stats)
            hsp.setAlignability(minValue)
                


