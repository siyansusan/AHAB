#class Target:
#    def __init__(self,ID,pos,seq):
#        self.ID=ID
#        self.pos=pos
#        self.seq=seq

#def loadTargets(filename):
#    targets=[]
#    with open(filename,"rt") as IN:
#        for line in IN:
#            fields=line.rstrip().split()
#            if(len(fields)!=3): continue
#            (ID,pos,seq)=fields
#            targets.append(Target(ID,int(pos),seq))
#    return targets


    #readGapLengths=anno.getReadGapLengths()
    #refGapLengths=anno.getRefGapLengths()
    #if(len(refGapLengths)>0 and max(refGapLengths)>MAX_REF_GAP): continue



#ALIGNABILITIES=open("alignabilities.txt","wt")
#REFGAPS=open("ref-gaps.txt","wt")



    #print(readGroup.getID(),"\t==> after clustering there are ",len(HSPs),
    #      " HSPs",sep="")
    #for hsp in HSPs:
    #    print("\t",hsp.toString(),sep="")
    #if(len(readGapLengths)>0):
    #    print("\tREAD GAPS:",",".join([str(x) for x in readGapLengths]))
    #if(len(refGapLengths)>0):
    #    print("\tREF GAPS:",",".join([str(x) for x in refGapLengths]))
    #    m=max(refGapLengths)
    #    print(m,file=REFGAPS,flush=True)
    #readsKept+=1
    #print("----------------------------------------------------------------")
#ALIGNABILITIES.close()
#REFGAPS.close()
#print(readsKept,"reads kept")

