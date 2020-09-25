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

def getAlignabilities(anno,bigwig):
    for hsp in anno.getHSPs():
        refCoords=hsp.getRefInterval()
        stats=bigwig.stats(hsp.getRefName(),refCoords.getBegin(),
                           refCoords.getEnd(),type="min")
        #for x in stats: print(x,file=OUT,flush=True)
        minValue=min(stats)
        hsp.setAlignability(minValue)

def process1HSP(anno):
    # Check whether this shows evidence of NO deletion, or whether it
    # is off-target
    pass

def process2HSPs(anno):
    HSPs=anno.getHSPs(); hsp1=HSPs[0]; hsp2=HSPs[1]
    interval1=hsp1.getRefInterval(); interval2=hsp2.getRefInterval()
    gaps=anno.getRefGaps()
    if(len(gaps)!=1): return
    gap=gaps[0]
    if(not gap.overlaps(DELETION_REGION)): 
        print("possible off-target edit:",anno.getReadID())
        return ### should report this as possible off-target edit

    d1=FIRST_CUT_SITE-interval1.getEnd()
    d2=interval2.getBegin()-SECOND_CUT_SITE
    #print("distances:",d1,d2,sep="\t")
    if(max(abs(d1),abs(d2))>MAX_ANCHOR_DISTANCE): 
        print("on-target but imperfect edit:",anno.getReadID())
        return ### should report this as on-target but imperfect edit
    print("on-target near-perfect deletion:",anno.getReadID())
    #MAX_ANCHOR_OVERLAP

def process3HSPs(anno):
    # Check for evidence of integration or inversion at target site, or of
    # off-target edits
    pass

def processManyHSPs(anno):
    # Complex edits: check whether on-target or off-target
    pass

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <filename.sam> <settings.config>\n")
(samFile,configFilename)=sys.argv[1:]

# Process configuration file
config=ConfigFile(configFilename)
dedup=config.lookupOrDie("DEDUPLICATE")
dedup=True if dedup=="True" or dedup=="true" or dedup=="T" or dedup=="TRUE" \
    or dedup=="Yes" or dedup=="yes" else False
MIN_IDENTITY=float(config.lookupOrDie("MIN_IDENTITY"))
MAX_REF_GAP=float(config.lookupOrDie("MAX_REF_GAP"))
alignabilityMap=config.lookupOrDie("ALIGNABILITY")
MIN_ALIGNABILITY=float(config.lookupOrDie("MIN_ALIGNABILITY"))

# Load target locations
#targets=loadTargets(config.lookupOrDie("TARGET_SITES"))
FIRST_CUT_SITE=int(config.lookupOrDie("FIRST_CUT_SITE"))
SECOND_CUT_SITE=int(config.lookupOrDie("SECOND_CUT_SITE"))
DELETION_REGION=Interval(FIRST_CUT_SITE,SECOND_CUT_SITE)
TARGET_CHROM=config.lookupOrDie("TARGET_CHROM")
MAX_ANCHOR_DISTANCE=int(config.lookupOrDie("MAX_ANCHOR_DISTANCE"))
MAX_ANCHOR_OVERLAP=int(config.lookupOrDie("MAX_ANCHOR_OVERLAP"))

# Open the ENCODE alignability map
bigwig=pyBigWig.open(alignabilityMap)

# Process SAM file
hspFactory=SamHspFactory()
stream=SamPairedReadStream(samFile)
#readsKept=0
#ALIGNABILITIES=open("alignabilities.txt","wt")
#REFGAPS=open("ref-gaps.txt","wt")
while(True):
    readGroup=stream.nextGroup()
    if(readGroup is None): break
    firstReads=getReadEnds(readGroup,1)

    # Should collapse these down to a single line based on a
    # SamAnnotationFactory:
    HSPs=hspFactory.makeHSPs(firstReads)
    HSPs=SamHspClusterer.cluster(HSPs)
    anno=SamAnnotation(HSPs)

    # Filter based on alignment quality and target chromosome
    if(not anno.allRefsSame()): continue
    if(anno.firstRef()!=TARGET_CHROM): continue
    if(anno.lowestPercentIdentity()<MIN_IDENTITY): continue
    getAlignabilities(anno,bigwig)
    if(anno.getLowestAlignability()<MIN_ALIGNABILITY): continue
    if(anno.anyRefsOverlap()): continue

    # Address cases of 1 HSP, 2 HSPs, 3 HSPs, and >3 HSPs
    numHSPs=anno.numHSPs()
    if(numHSPs==1): process1HSP(anno)
    elif(numHSPs==2): process2HSPs(anno)
    elif(numHSPs==3): process3HSPs(anno)
    else: processManyHSPs(anno)

    #if(not anno.allSameStrand()): continue

    #readGapLengths=anno.getReadGapLengths()
    #refGapLengths=anno.getRefGapLengths()
    #if(len(refGapLengths)>0 and max(refGapLengths)>MAX_REF_GAP): continue
    

