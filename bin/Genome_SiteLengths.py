#!/usr/bin/env python3

from Bio import SeqIO
import os
import sys

path=sys.argv[1]
contigFile=(path+'/contigs.fa')

file = open(path+'/contigs_SeqLength.tsv', "w") 

for seq_record in SeqIO.parse(contigFile,"fasta"):
    file.write(str(seq_record.id)+"\t"+str(len(seq_record))+"\n")
file.close()

#Create file with an entry for each individual site in the alignment
siteCount=0
locListFile = open(path+'/contigs_LocList','a+')
with open(path +"/contigs_SeqLength.tsv","r") as filein:
    for line in iter(filein):
        splitline=line.split()
        for x in range(1,(int(splitline[1])+1)):
            locListFile.write((splitline[0] +'/'+str(x)+'\n'))
            siteCount+=1
locListFile.close()

