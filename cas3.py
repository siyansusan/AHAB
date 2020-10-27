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
#  * Let user specify -1 as "N/A" for any parameter
#  * We're only processing the first read of each pair
#  * Add comments here and in all support classes
#=========================================================================


#=========================================================================
#                              class Cas3Analysis
#=========================================================================
class Cas3Analysis(Ahab):
    def __init__(self,configFile,outputDir):
        super().__init__(configFile)
        config=self.config
        dedup=config.lookupOrDie("DEDUPLICATE")
        dedup=True if dedup in ("True","true","T","TRUE","Yes","yes") \
            else False
        self.dedup=dedup
        self.MIN_IDENTITY=float(config.lookupOrDie("MIN_IDENTITY"))
        self.MAX_REF_GAP=float(config.lookupOrDie("MAX_REF_GAP"))
        self.MIN_ALIGNABILITY=float(config.lookupOrDie("MIN_ALIGNABILITY"))
        self.CUT_SITE=int(config.lookupOrDie("CUT_SITE"))
        self.TARGET_CHROM=config.lookupOrDie("TARGET_CHROM")
        self.PRIMER_SEQ=config.lookupOrDie("PRIMER_SEQ")
        #self.MAX_ANCHOR_DISTANCE=int(config.lookupOrDie("MAX_ANCHOR_DISTANCE"))
        #self.MAX_ANCHOR_OVERLAP=int(config.lookupOrDie("MAX_ANCHOR_OVERLAP"))
        self.MAX_READ_GAP=int(config.lookupOrDie("MAX_READ_GAP"))
        self.MIN_ALIGNED_PROPORTION=\
            float(config.lookupOrDie("MIN_ALIGNED_PROPORTION"))
        self.OUTPUT_DIR=outputDir
        self.prepareOutputFiles(self.OUTPUT_DIR)
        self.CHROMS=set(("chr1","chr2","chr3","chr4","chr5","chr6",
                         "chr7","chr8","chr9","chr10","chr11","chr12",
                         "chr13","chr14","chr15","chr16","chr17","chr18",
                         "chr19","chr20","chr21","chr22","chrX","chrY"))

    def __del__(self):
        self.ON_TARGET_DELETION.close()
        self.ON_TARGET_IMPERFECT_DELETION.close()
        self.ON_TARGET_NO_EDIT.close()
        self.ON_TARGET_SHORT_INDELS.close()
        self.OFF_TARGET_EDIT.close()
        self.OFF_TARGET_NO_EDIT.close()
        self.MISC.close()
        self.FAILED_FILTER.close()

    def prepareOutputFiles(self,DIR):
        if(not os.path.exists(DIR)):
            os.system("mkdir -p "+DIR)
        if(not rex.find("(.*)/$",DIR)): DIR=DIR+"/"
        self.ON_TARGET_DELETION=open(DIR+"bin-on-target-deletion.txt","wt")
        self.ON_TARGET_IMPERFECT_DELETION=\
               open(DIR+"bin-on-target-imperfect-deletion.txt","wt")
        self.ON_TARGET_NO_EDIT=open(DIR+"bin-on-target-no-edit.txt","wt")
        self.ON_TARGET_SHORT_INDELS=\
            open(DIR+"bin-on-target-short-indels","wt")
        self.OFF_TARGET_EDIT=open(DIR+"bin-off-target-edit.txt","wt")
        self.OFF_TARGET_NO_EDIT=open(DIR+"bin-off-target-no-edit.txt","wt")
        self.MISC=open(DIR+"bin-misc.txt","wt")
        self.FAILED_FILTER=open(DIR+"bin-failed-filter.txt","wt")

    def processCases(self,anno):
        numHSPs=anno.numHSPs()
        if(numHSPs==1): self.process1HSP(anno)
        elif(numHSPs==2): self.process2HSPs(anno)
        elif(numHSPs==3): self.process3HSPs(anno)
        else: self.processManyHSPs(anno)

    def process1HSP(self,anno):
        # Precondition: 1 HSP only, and is on the target chromosome

        hsp=anno.getHSPs()[0]
        if(self.PRIMER_SEQ in anno.getSamRecord().getSequence()):
            if(hsp.containsIndels()):
                self.bin(anno,self.ON_TARGET_SHORT_INDELS)
            else:
                self.bin(anno,self.ON_TARGET_NO_EDIT)
            return
        if(hsp.containsIndels()):
            self.bin(anno,self.OFF_TARGET_EDIT)
        else:
            self.bin(anno,self.OFF_TARGET_NO_EDIT)


    def process1HSP_OLD(self,anno):
        # Precondition: 1 HSP only, and is on the target chromosome
        hsp=anno.getHSPs()[0]
        interval=hsp.getRefInterval()
        if(interval.contains(self.CUT_SITE)):
            if(hsp.containsIndels()):
                self.bin(anno,self.ON_TARGET_SHORT_INDELS)
            else:
                self.bin(anno,self.ON_TARGET_NO_EDIT)
            return
        if(hsp.containsIndels()):
            self.bin(anno,self.OFF_TARGET_EDIT)
        else:
            self.bin(anno,self.OFF_TARGET_NO_EDIT)

    def process2HSPs(self,anno):
        if(anno.allSameStrand()): 
            self.process2HSPsSameStrand(anno)
        else:
            self.process2HSPsDifferentStrand(anno)

    def process2HSPsDifferentStrand(self,anno):
        self.bin(anno,self.MISC)

    def checkRefGapDeletion(self,anno):
        gaps=anno.getRefGaps()
        if(len(gaps)!=1): return False
        gap=gaps[0]
        if(gap.getLength()>self.MAX_REF_GAP):
            print("XXX1",gap.getLength(),"> MAX_REF_GAP") ###
            self.bin(anno,self.FAILED_FILTER)
            return False
        if(self.PRIMER_SEQ not in anno.getSamRecord().getSequence()):
            print("XXX10 NO PRIMER SEQ")
        #if(not gap.contains(self.CUT_SITE)): 
        #    print("XXX2",gap.toString(),"does not contain",self.CUT_SITE) ###
            self.bin(anno,self.OFF_TARGET_EDIT)
            return False
        return True

    def readGapSmallerThan(self,anno,MAX):
        gaps=anno.getReadGapLengths()
        if(len(gaps)==0): return True
        if(gaps[0]>=MAX): print("XXX3",gaps[0],">",MAX) ###
        return gaps[0]<MAX
            
    def process2HSPsSameStrand(self,anno):
        HSPs=anno.getHSPs(); hsp1=HSPs[0]; hsp2=HSPs[1]
        interval1=hsp1.getRefInterval(); interval2=hsp2.getRefInterval()
        if(not self.checkRefGapDeletion(anno)): return

        # Check some more filters
        if(not self.readGapSmallerThan(anno,self.MAX_READ_GAP) or
           anno.alignedProportion()<self.MIN_ALIGNED_PROPORTION):
            if(anno.alignedProportion()<self.MIN_ALIGNED_PROPORTION):
                print("XXX4",anno.alignedProportion(),"<",self.MIN_ALIGNED_PROPORTION)
            self.bin(anno,self.FAILED_FILTER)
            return
        
        self.bin(anno,self.ON_TARGET_DELETION)

    def process3HSPs(self,anno):
        ### Need to check for evidence of integration or inversion at target 
        ### site, or of off-target edits
        self.bin(anno,self.MISC)
        pass

    def processManyHSPs(self,anno):
        ### Complex edits: need to check whether on-target or off-target
        self.bin(anno,self.MISC)

    def filter(self,anno):
        if(not anno.allRefsSame()): 
            print("XXX5")
            self.bin(anno,self.FAILED_FILTER)
            return False
        if(anno.firstRef()!=ahab.TARGET_CHROM): 
            print("XXX6",anno.firstRef(),"!=",ahab.TARGET_CHROM)
            self.bin(anno,self.FAILED_FILTER)
            return False
        if(anno.lowestPercentIdentity()<ahab.MIN_IDENTITY): 
            print("XXX7",anno.lowestPercentIdentity(),"<",ahab.MIN_IDENTITY)
            self.bin(anno,self.FAILED_FILTER)
            return False
        ahab.getAlignabilities(anno)
        if(anno.getLowestAlignability()<ahab.MIN_ALIGNABILITY):
            print("XXX8",anno.getLowestAlignability(),"<",ahab.MIN_ALIGNABILITY)
            self.bin(anno,self.FAILED_FILTER)
            return False
        if(anno.anyRefsOverlap()): 
            print("XXX9")
            self.bin(anno,self.FAILED_FILTER)
            return False
        return True

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <settings.config> <filename.sam> <output-dir>\n")
(configFilename,samFile,outputDir)=sys.argv[1:]

# Instantiate Ahab object
ahab=Cas3Analysis(configFilename,outputDir)

# Process SAM file
hspFactory=SamHspFactory()
stream=SamPairedReadStream(samFile)
readsSeen=0
while(True):
    readGroup=stream.nextGroup()
    if(readGroup is None): break
    readsSeen+=1
    firstReads=readGroup.getReadEnds(1)

    ### NOTE: we're currently only processing the first read of each pair!

    # Cluster the HSPs and produce an Alignment object
    HSPs=hspFactory.makeHSPs(firstReads)
    HSPs=SamHspClusterer.cluster(HSPs)
    anno=SamAnnotation(HSPs)

    # Filter based on alignment quality and target chromosome, etc.
    if(not ahab.filter(anno)): continue

    # Address cases of 1 HSP, 2 HSPs, 3 HSPs, and >3 HSPs
    ahab.processCases(anno)

print(ahab.readsBinned,"reads binned, out of ",readsSeen)    
del ahab # Call destructor to clean up

