# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:54:20 2020

@author: mcveigh
"""

#Feature Slicer to extract mito small subunit rRNA sequences from a larger file

import Bio
import sys
#from Bio import SeqIO, SeqFeature
import os

inputfile = sys.argv[1]
outputfile = sys.argv[2]

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
sequence = []  
for seq_record in SeqIO.parse(inputfile, "genbank"):    
    str_id = seq_record.id
    counter = 0
    for feature in seq_record.features:
        if feature.type == 'rRNA':
            #counter = 0
            seq_record.description = seq_record.annotations["organism"]
            rnaname = str(feature.qualifiers.get("product"))   
            if rnaname == "['16S ribosomal RNA']" or rnaname == "['small subunit ribosomal RNA']": 
                counter += 1
                if counter == 1:
                    #mystart = feature.location._start.position + 1
                    #myend = feature.location._end.position
                    #print("I found a 23S rRNA in", seq_record.id, "from", mystart, "to", myend)
                    if feature.strand == -1:
                        sub_record = seq_record[feature.location.start:feature.location.end].reverse_complement()
                        sub_record.id = str_id
                        sub_record.description = (seq_record.description + " small subunit ribosomal RNA")
                    elif feature.strand == 1:
                        sub_record = seq_record[feature.location.start:feature.location.end]
                        sub_record.description = (seq_record.description + " small subunit ribosomal RNA")
                        print(sub_record)
                        
                #else:
                if counter == 1:
                    sequence.append(sub_record)
                    print(sub_record.id)
                elif counter > 1:
                    print(str(str_id)+" has "+str(counter)+" small subunit ribosomal RNA features")
                else:
                    print(str(str_id)+" has no small subunit ribosomal RNA features")

SeqIOwrite(sequence, outputfile, "fasta")  .
    
fh = open(outputfile)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
fh.close()
print("The number of sequences in the final output file is: ", n)
print ("FeatureSlicer script is done and output is saved in:", outputfile)   
    
    

    