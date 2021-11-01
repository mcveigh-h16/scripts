# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:32:45 2020

@author: mcveigh
"""
#Feature Slicer to extract 23S rRNA sequences from a larger file

import Bio
import sys
#from Bio import SeqIO, SeqFeature
import os

inputfile = sys.argv[1]
outputfile = sys.argv[2]

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
sequence = []  
for seq_record in SeqIO.parse(inputfile, "genbank"):    
    str_id = seq_record.id
    counter = 0
    for feature in seq_record.features:
        if feature.type == 'CDS':
            #counter = 0
            seq_record.description = seq_record.annotations["organism"]
            proteinname = str(feature.qualifiers.get("product"))   
            if proteinname == "['cytochrome b']":
                counter += 1
                if counter == 1:
                    mystart = feature.location._start.position + 1
                    myend = feature.location._end.position
                    print("I found a cytb in", seq_record.id, "from", mystart, "to", myend)
                    if feature.strand == -1:
                        sub_record = seq_record[feature.location.start:feature.location.end].reverse_complement()
                        sub_record.id = str_id
                        sub_record.description = (seq_record.description + " cytochrome b")
                    elif feature.strand == 1:
                            sub_record = seq_record[feature.location.start:feature.location.end]
                            sub_record.description = (seq_record.description + " cytochrome b")
                #else:
            if counter == 1:
                sequence.append(sub_record)
            elif counter > 1:
                print(str(str_id)+" has "+str(counter)+" cytochrome b features")
            elif counter == 0:
                print(str(str_id)+" has no cytochrome b features")
        else:

SeqIO.write(sequence, outputfile, "fasta")  
    
fh = open(outputfile)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
fh.close()
print("The number of sequences in the final output file is: ", n)
print ("FeatureSlicer script is done and output is saved in:", outputfile)   
    
    

    