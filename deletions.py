#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import ProgramName
from SamReader import SamReader
from FastaReader import FastaReader
from FastaWriter import FastaWriter
from Translation import Translation
from Interval import Interval
from Pipe import Pipe
from Rex import Rex
rex=Rex()
import TempFilename
from CigarString import CigarString
from SamPairedReadStream import SamPairedReadStream
from SamHspFactory import SamHspFactory
from SamHspClusterer import SamHspClusterer
from SamAnnotation import SamAnnotation

MIN_IDENTITY=0.9

class Target:
    def __init__(self,ID,pos,seq):
        self.ID=ID
        self.pos=pos
        self.seq=seq

def loadTargets(filename):
    targets=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=3): continue
            (ID,pos,seq)=fields
            targets.append(Target(ID,int(pos),seq))
    return targets

def getReadEnds(readGroup,end1or2):
    reads=[]
    for pair in readGroup.getReads():
        if(end1or2==1): reads.append(pair.read1)
        else: reads.append(pair.read2)
    return reads

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <filename.sam> <target-sites.txt> <deduplicate:0/1>\n")
(samFile,targetFile,dedup)=sys.argv[1:]
dedup=True if dedup==1 else False

# Load target locations
targets=loadTargets(targetFile)

# Process SAM file
hspFactory=SamHspFactory()
stream=SamPairedReadStream(samFile)
while(True):
    readGroup=stream.nextGroup()
    if(readGroup is None): break
    firstReads=getReadEnds(readGroup,1)

    # Should collapse these down to a single line based on a
    # SamAnnotationFactory:
    HSPs=hspFactory.makeHSPs(firstReads)
    HSPs=SamHspClusterer.cluster(HSPs)
    anno=SamAnnotation(HSPs)

    numHSPs=anno.numHSPs()
    if(not anno.allRefsSame()): continue ### debugging
    if(anno.firstRef()!="chrX"): continue ### debugging
    if(numHSPs<2): continue
    if(not anno.allSameStrand()): continue
    if(anno.lowestPercentIdentity()<MIN_IDENTITY): continue

    print(readGroup.getID(),"\t==> after clustering there are ",len(HSPs),
          " HSPs",sep="")
    for hsp in HSPs:
        print(hsp.toString(),"\t",sep="",end="")
    print()

    readGapLengths=anno.getReadGapLengths()
    refGapLengths=anno.getRefGapLengths()
    if(len(readGapLengths)>0):
        print("\tREAD GAPS:",",".join([str(x) for x in readGapLengths]))
    if(len(refGapLengths)>0):
        print("\tREF GAPS:",",".join([str(x) for x in refGapLengths]))
        if(max(refGapLengths)>1000000): 
            print("\tMAX GAP LENGTH EXCEEDED")

    
