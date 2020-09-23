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
import pyBigWig
from ConfigFile import ConfigFile

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

def getAlignabilities(anno,bigwig,OUT):
    for hsp in anno.getHSPs():
        refCoords=hsp.getRefInterval()
        stats=bigwig.stats(hsp.getRefName(),refCoords.getBegin(),
                           refCoords.getEnd(),type="min")
        for x in stats: print(x,file=OUT,flush=True)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <filename.sam> <settings.config>\n")
(samFile,configFilename)=sys.argv[1:]

# Process configuration file
config=ConfigFile(configFilename)
dedup=configFile.lookupOrDie("DEDUPLICATE")
dedup=True if dedup=="True" or dedup=="true" or dedup=="T" or dedup=="TRUE" \
    or dedup=="Yes" or dedup=="yes" else False
MIN_IDENTITY=float(configFile.lookupOrDie("MIN_IDENTITY"))
MAX_REF_GAP=float(configFile.lookupOrDie("MAX_REF_GAP"))
alignabilityMap=configFile.lookupOrDie("ALIGNABILITY")

# Load target locations
targets=loadTargets(configFile.lookupOrDie("TARGET_SITES"))

# Open the ENCODE alignability map
bigwig=pyBigWig.open(alignabilityMap)

# Process SAM file
hspFactory=SamHspFactory()
stream=SamPairedReadStream(samFile)
#numErrors=0
numGreater1000000=0; numGreater100000=0; numGreater10000=0; numGreater1000=0
readsKept=0
ALIGNABILITIES=open("alignabilities.txt","wt")
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
    if(anno.anyRefsOverlap()): continue
    if(not anno.allSameStrand()): continue
    if(anno.lowestPercentIdentity()<MIN_IDENTITY): continue
    getAlignabilities(anno,bigwig,ALIGNABILITIES)

    readGapLengths=anno.getReadGapLengths()
    refGapLengths=anno.getRefGapLengths()
    #if(len(refGapLengths)==0 or max(refGapLengths)<MAX_REF_GAP): continue ###

    print(readGroup.getID(),"\t==> after clustering there are ",len(HSPs),
          " HSPs",sep="")
    for hsp in HSPs:
        print("\t",hsp.toString(),sep="")
    if(len(readGapLengths)>0):
        print("\tREAD GAPS:",",".join([str(x) for x in readGapLengths]))
    if(len(refGapLengths)>0):
        print("\tREF GAPS:",",".join([str(x) for x in refGapLengths]))
        m=max(refGapLengths)
        if(m>1000000): numGreater1000000+=1
        elif(m>100000): numGreater100000+=1
        elif(m>10000): numGreater10000+=1
        elif(m>1000): numGreater1000+=1
        #if(max(refGapLengths)>MAX_REF_GAP): 
        #    print("\tMAX GAP LENGTH EXCEEDED")
            #numErrors+=1
    readsKept+=1
    print("----------------------------------------------------------------")
ALIGNABILITES.close()

#print(numErrors,"total errors")
print(numGreater1000000,">Mb")
print(numGreater100000,">100k")
print(numGreater10000,">10k")
print(numGreater1000,">1k")
print("out of",readsKept,"reads kept")
