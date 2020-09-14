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
    for pair in readGroup:
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
#reader=SamReader(samFile)
hspFactory=SamHspFactory()
stream=SamPairedReadStream(samFile)
while(True):
    readGroup=stream.nextGroup()
    if(readGroup is None): break
    firstReads=getReadEnds(readGroup,1)
    HSPs=hspFactory.makeHSPs(firstReads)
    print(len(readGroup),"reads in group",len(HSPs),"HSPs")
    HSPs=SamHspClusterer.cluster(HSPs)
    print("\t==> after clustering there are ",len(HSPs),"HSPs")

#while(True):
#    pair=stream.nextPair()
#    if(pair is None): break
#    read1=pair.read1; read2=pair.read2
#    score=pair.computeScore()
#    print(pair.getID(),score,read1.getCigar().toString(),
#          read2.getCigar().toString(),"".join(read1.parseMDtag()),
#          "".join(read1.parseMDtag()),sep="\t")

    #print("pair:",read1.getID(),read2.getID(),sep="\t")
    #if(read1.getRefName()!=read2.getRefName()):
    #    print("REF MISMATCH IN PAIRED READ:",read1.getRefName(),
    #          read2.getRefName(),sep="\t")
    #continue
    #rec=reader.nextSequence()
    #if(rec is None): break
    #if(rec.flag_unmapped()): continue
    #if(dedup and rec.flag_PCRduplicate()): continue
    #firstOfPair=rec.flag_firstOfPair()
    #cigar=rec.CIGAR
    #if(cigar.completeMatch()): continue
    #cigar.computeIntervals(rec.getRefPos())
    #readLen=len(rec.getSequence())
    #MDfields=rec.parseMDtag()
    #refName=rec.getRefName()
    #refPos=rec.getRefPos()
    #print(rec.ID,firstOfPair,refName,refPos,cigar.toString(),MDfields,sep="\t")



