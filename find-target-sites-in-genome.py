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
import ProgramName
from FastaReader import FastaReader
from Translation import Translation

V_50_9="CACCACTCACCTC"
V_51_1="GACCATTTCCCAC"
V_50_9=Translation.reverseComplement(V_50_9)
V_51_1=Translation.reverseComplement(V_51_1)

CHR_X="/home/bmajoros/Reference_Data/Genomes/hg19/chrX.fa"

(defline,chrX)=FastaReader.firstSequence(CHR_X)
for i in range(31140035,33229429):
    if(V_50_9==chrX[i:(i+13)]): print("V_50_9\t",i)
    if(V_51_1==chrX[i:(i+13)]): print("V_51_1\t",i+13)







