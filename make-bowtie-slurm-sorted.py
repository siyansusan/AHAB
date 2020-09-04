#!/usr/bin/python3

import ProgramName
import sys
import os
import subprocess
import re
from Rex import Rex
from SlurmWriter import SlurmWriter


rex=Rex()
writer = SlurmWriter()

if(len(sys.argv)!=5):
    exit(ProgramName.get()+ " <fastq-dir> <alignment-dir> <working-dir> <reference>\n")
(fastqDir,alignmentDir,workingDir,ref)=sys.argv[1:]

if(rex.find("\/$",workingDir)): workingDir=workingDir[:-1]
if(rex.find("\/$",fastqDir)): fastqDir=fastqDir[:-1]
#match=re.search("\/$",workingDir)
slurmDir=workingDir+"/bowtie-slurms"


isDir = os.path.isdir(slurmDir)
if isDir:
    proc = subprocess.Popen("rm -f -r " + slurmDir + "/*", shell=True, executable="/bin/bash")
    proc.communicate()

else:
    proc = subprocess.Popen("mkdir "+str(slurmDir), shell=True, executable="/bin/bash")
    proc.communicate()

# not sure where the directory is
#path=fastqDir/*"FWD.fastq.gz"
files = os.listdir(fastqDir) #fastq $ for pl



for file in files:
    if(rex.find("([^/]+)\.FWD\.fq\.gz$",file.rstrip())):
        libraryID=rex[1]
        file1= fastqDir+"/"+libraryID+".FWD.fq.gz"
        file2= fastqDir+"/"+libraryID+".REV.fq.gz"
        print(file1)
        if not os.path.isfile(file1):
            raise Exception("FWD file not found")
        if not os.path.isfile(file2):
            raise Exception("REV file not found")
        #CPU 32
        cmd = [
            f"cd {workingDir}",
            f"module load bowtie2",
            f"bowtie2 --local -k 20 -p 32 -x {ref} -1 {file1} -2 {file2} | samtools sort -@32 -O SAM -n > {alignmentDir}/{libraryID}.sam"
        ]
        cmd = " ; ".join(cmd)
        #cmd = "cd "+workingDir+" ; module load bowtie2 ; bowtie2 --local -k 20 -p 32 -x ../../bill/bowtie/exon51 -1 "+file1+" -2 "+file2+" | samtools view -@ 32 -b > "+alignmentDir+"/"+libraryID+".bam"
        writer.addCommand(cmd)


writer.mem(32000)
writer.setQueue("all")
writer.writeArrayScript(
    slurmDir,
    "BOWTIE",
    300,
    ""
    )






