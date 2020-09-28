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
from Interval import Interval
from Rex import Rex
rex=Rex()
from SamPairedReadStream import SamPairedReadStream
from SamHspFactory import SamHspFactory
from SamHspClusterer import SamHspClusterer
from SamAnnotation import SamAnnotation
from Ahab import Ahab

#=========================================================================
# TODO LIST:
#  * We're only processing the first read of each pair
#  * 
#=========================================================================


#=========================================================================
#                              class Analysis
#=========================================================================
class Analysis(Ahab):
    def __init__(self,configFile):
        super().__init__(configFile)
        config=self.config
        dedup=config.lookupOrDie("DEDUPLICATE")
        dedup=True if dedup in ("True","true","T","TRUE","Yes","yes") \
            else False
        self.dedup=dedup
        self.MIN_IDENTITY=float(config.lookupOrDie("MIN_IDENTITY"))
        self.MAX_REF_GAP=float(config.lookupOrDie("MAX_REF_GAP"))
        self.MIN_ALIGNABILITY=float(config.lookupOrDie("MIN_ALIGNABILITY"))
        self.FIRST_CUT_SITE=int(config.lookupOrDie("FIRST_CUT_SITE"))
        self.SECOND_CUT_SITE=int(config.lookupOrDie("SECOND_CUT_SITE"))
        self.DELETION_REGION=Interval(self.FIRST_CUT_SITE,self.SECOND_CUT_SITE)
        self.TARGET_CHROM=config.lookupOrDie("TARGET_CHROM")
        self.MAX_ANCHOR_DISTANCE=int(config.lookupOrDie("MAX_ANCHOR_DISTANCE"))
        #self.MAX_ANCHOR_OVERLAP=int(config.lookupOrDie("MAX_ANCHOR_OVERLAP"))
        self.OUTPUT_DIR=config.lookupOrDie("OUTPUT_DIR")
        self.prepareOutputFiles(self.OUTPUT_DIR)

    def __del__(self):
        self.ON_TARGET_DELETION.close()
        self.ON_TARGET_IMPERFECT_DELETION.close()
        self.ON_TARGET_NO_EDIT.close()
        self.ON_TARGET_INVERSION.close()
        self.OFF_TARGET_EDIT.close()
        self.OFF_TARGET_NO_EDIT.close()

    def prepareOutputFiles(self,DIR):
        if(not os.path.exists(DIR)):
            os.system("mkdir -p "+DIR)
        if(not rex.find("(.*)/$",DIR)): DIR=DIR+"/"
        self.ON_TARGET_DELETION=open(DIR+"bin-on-target-deletion.txt","wt")
        self.ON_TARGET_IMPERFECT_DELETION=\
               open(DIR+"bin-on-target-imperfect-deletion.txt","wt")
        self.ON_TARGET_NO_EDIT=open(DIR+"bin-on-target-no-edit.txt","wt")
        self.ON_TARGET_INVERSION=open(DIR+"bin-on-target-inversion.txt","wt")
        self.OFF_TARGET_EDIT=open(DIR+"bin-off-target-edit.txt","wt")
        self.OFF_TARGET_NO_EDIT=open(DIR+"bin-off-target-no-edit.txt","wt")

    def processCases(self,anno):
        numHSPs=anno.numHSPs()
        if(numHSPs==1): self.process1HSP(anno)
        elif(numHSPs==2): self.process2HSPs(anno)
        elif(numHSPs==3): self.process3HSPs(anno)
        else: self.processManyHSPs(anno)

    def process1HSP(self,anno):
        # Check whether this shows evidence of NO deletion, or whether it
        # is off-target
        pass

    def process2HSPs(self,anno):
        if(anno.allSameStrand()): 
            self.process2HSPsSameStrand(anno)
        else:
            self.process2HSPsDifferentStrand(anno)

    def process2HSPsDifferentStrand(self,anno):
        # Check for inversions
        HSPs=anno.getHSPs(); hsp1=HSPs[0]; hsp2=HSPs[1]
        interval1=hsp1.getRefInterval(); interval2=hsp2.getRefInterval()
        FIRST_CUT_SITE=self.FIRST_CUT_SITE
        SECOND_CUT_SITE=self.SECOND_CUT_SITE
        MAX_ANCHOR_DISTANCE=self.MAX_ANCHOR_DISTANCE
        if(abs(FIRST_CUT_SITE-interval1.getEnd())<=MAX_ANCHOR_DISTANCE and
           abs(SECOND_CUT_SITE-interval2.getEnd())<=MAX_ANCHOR_DISTANCE or
           abs(interval1.getBegin()-FIRST_CUT_SITE)<=MAX_ANCHOR_DISTANCE and
           abs(interval2.getBegin()-SECOND_CUT_SITE)<=MAX_ANCHOR_DISTANCE):
            self.bin(anno.getReadID(),self.ON_TARGET_INVERSION)
            
    def process2HSPsSameStrand(self,anno):
        HSPs=anno.getHSPs(); hsp1=HSPs[0]; hsp2=HSPs[1]
        interval1=hsp1.getRefInterval(); interval2=hsp2.getRefInterval()
        gaps=anno.getRefGaps()
        if(len(gaps)!=1): return
        gap=gaps[0]
        if(not gap.overlaps(self.DELETION_REGION)): 
            self.bin(anno.getReadID(),self.OFF_TARGET_EDIT)
            return
        FIRST_CUT_SITE=self.FIRST_CUT_SITE
        SECOND_CUT_SITE=self.SECOND_CUT_SITE
        MAX_ANCHOR_DISTANCE=self.MAX_ANCHOR_DISTANCE
        d1=FIRST_CUT_SITE-interval1.getEnd()
        d2=interval2.getBegin()-SECOND_CUT_SITE
        if(max(abs(d1),abs(d2))>MAX_ANCHOR_DISTANCE): 
            self.bin(anno.getReadID(),self.ON_TARGET_IMPERFECT_DELETION)
            return
        self.bin(anno.getReadID(),self.ON_TARGET_DELETION)

    def process3HSPs(self,anno):
        # Check for evidence of integration or inversion at target site, or of
        # off-target edits
        pass

    def processManyHSPs(self,anno):
        # Complex edits: check whether on-target or off-target
        pass

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <filename.sam> <settings.config>\n")
(samFile,configFilename)=sys.argv[1:]

# Instantiate Ahab object
ahab=Analysis(configFilename)

# Process SAM file
hspFactory=SamHspFactory()
stream=SamPairedReadStream(samFile)
readsKept=0
while(True):
    readGroup=stream.nextGroup()
    if(readGroup is None): break
    firstReads=readGroup.getReadEnds(1)

    ### NOTE: we're currently only processing the first read of each pair!

    # Should collapse these down to a single line based on a
    # SamAnnotationFactory:
    HSPs=hspFactory.makeHSPs(firstReads)
    HSPs=SamHspClusterer.cluster(HSPs)
    anno=SamAnnotation(HSPs)

    # Filter based on alignment quality and target chromosome
    if(not anno.allRefsSame()): continue
    if(anno.firstRef()!=ahab.TARGET_CHROM): continue
    if(anno.lowestPercentIdentity()<ahab.MIN_IDENTITY): continue
    ahab.getAlignabilities(anno)
    if(anno.getLowestAlignability()<ahab.MIN_ALIGNABILITY): continue
    if(anno.anyRefsOverlap()): continue

    # Address cases of 1 HSP, 2 HSPs, 3 HSPs, and >3 HSPs
    ahab.processCases(anno)

print(readsKept,"reads kept")    
del ahab # Call destructor to clean up

