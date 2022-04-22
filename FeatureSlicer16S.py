# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 13:32:45 2022

@author: mcveigh
"""
#Feature Slicer to download genomes from an accession list extract 16S rRNA sequences and determine if the 16S genes in a single genome are identical. 
#Identical sequences are disgarded from final outout. 

import Bio
import sys
#from Bio import SeqIO, SeqFeature
import os

inputfile = sys.argv[1]
outputfile = sys.argv[2]

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#os.system("efetch -db nuccore -input acclist -format gb > genome.gbk")

sequence = []  
for seq_record in SeqIO.parse(inputfile, "genbank"):    
    str_id = seq_record.id
    counter = 0
    for feature in seq_record.features:
        if feature.type == 'rRNA':
            #counter = 0
            seq_record.description = seq_record.annotations["organism"]
            rnaname = str(feature.qualifiers.get("product"))   
            if rnaname == "['16S ribosomal RNA']":
                counter += 1
                start = []
                stop = []
                mystart = feature.location._start.position + 1
                myend = feature.location._end.position
                
                if feature.strand == -1:
                    sub_record = seq_record[feature.location.start:feature.location.end].reverse_complement()
                    sub_record.id = str_id
                    print(seq_record.id, "minus strand", myend, "to", mystart)
                    start = str(myend)
                    stop = str(mystart)
                    sub_record.description = (seq_record.description + " " + start + " " + stop)
                    #sub_record.description = (seq_record.description + " 16S ribosomal RNA")
                elif feature.strand == 1:
                    start = str(mystart)
                    stop = str(myend)
                    sub_record = seq_record[feature.location.start:feature.location.end]
                    print(seq_record.id, "plus strand", mystart, "to", myend)
                    sub_record.description = (seq_record.description + " " + start + " " + stop)
                    #sub_record.description = (seq_record.description + " 16S ribosomal RNA")
                sequence.append(sub_record)
    if counter > 1:
        print(str(str_id)+" has "+str(counter)+" 16S ribosomal RNA features")
    else:
        print(str(str_id)+" has no 16S ribosomal RNA features")


SeqIO.write(sequence, outputfile, "fasta")  
    
fh = open(outputfile)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
fh.close()
print("The number of sequences in the final output file is: ", n)
print ("FeatureSlicer script is done and output is saved in:", outputfile)   
    
    

    