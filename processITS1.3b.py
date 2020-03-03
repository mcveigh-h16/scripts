# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 08:39:05 2020

@author: mcveigh
"""

#
# processITS designed to validate ITS sequences. Sequences that pass all tests are saved to a outfile file in fasta format
#Selectively sort GenBank flatfiles removing selected seqeunces and saving them to a new fasta file
# User must specify the input filename and the outputfile name

import pandas as pd
import Bio
import os
import sys

inputfile = sys.argv[1]
outputfile = sys.argv[2]
#print(inputfile)
#print(outputfile)

#Read in the reject list and find the accession
reject_file_name = (r'ITS_reject_seqs3.txt') 
reject_file_name = (r'/panfs/pan1.be-md.ncbi.nlm.nih.gov/dnaorg/ITS/ITS_reject_seqs')
rejectlist_df = pd.read_csv(reject_file_name, sep='\t', index_col=None, low_memory=False, header=None, names=["accession", "type", "reason"])
rejectlist = rejectlist_df['accession']
print (rejectlist) 
   
#Parse the GenBank file and remove any sequences found on the reject list
from Bio import SeqIO
sequences = [] 
found = []
#missingRNA = []

for seq_record in SeqIO.parse(inputfile, "genbank"):    
    str_id = seq_record.id
    #print(seq_record.id)
    if str_id.find('.') != -1:        
        str_id = str_id[:str_id.find('.')]        
    if str_id not in rejectlist_df['accession'].tolist():                       
        seq_record.description = seq_record.annotations["organism"]
        sequences.append(seq_record)
    else:
        found.append(seq_record)
        print("I found this accession on the reject list and wrote the sequence to found.fsa: ", seq_record.id)  
    SeqIO.write(sequences, "stripped.fsa", "fasta")  
    SeqIO.write(found, "found.fsa", "fasta")   
  
#Run short version of ribodbmaker.pl    
#os.system("ribodbmaker.pl -f --skipfribo1 --skipfribo2 --skipfmspan --skipingrup --skipclustr --skiplistms --skipmstbl --taxin /panfs/pan1/dnaorg/rrna/git-ncbi-rrna-project/taxonomy-files/ncbi-taxonomy-tree.ribodbmaker.txt stripped.fsa ribo-out")

#Run CMscan using the output of ribodbmaker.pl
#os.system("cmscan --cpu 4 --mid -T 20 --verbose --tblout tblout.df.txt rrna.cm ribo-out/ribo-out.ribodbmaker.final.fa > cmscanOUTPUT.df.txt")
#os.system("cmscan --cpu 4 --mid -T 20 --verbose --anytrunc --tblout tblout.at.txt rrna.cm ribo-out/ribo-out.ribodbmaker.final.fa > cmscanOUTPUT.at.txt")
os.system("cat tblout.df.txt tblout.at.txt > tblout.both.txt")
os.system("perl cmsearch_tblout_deoverlap/cmsearch-deoverlap.pl --maxkeep -s --cmscan tblout.both.txt")
os.system("head -n2 tblout.both.txt > final.tblout")
os.system("cat tblout.both.txt.deoverlapped >> final.tblout")

#Add seq length to cmscan output
os.system("esl-seqstat -a ribo-out/ribo-out.ribodbmaker.final.fa | grep ^\\= | awk '{ printf(\"%s %s\\n\", $2, $3); }' > my.seqlen")
os.system("perl tblout-add.pl -t final.tblout 18 my.seqlen 3 > cmscan_final.tblout")

#Parse the final results of CMscan, sort and write fasta files
CMscan_output = (r'cmscan_final.tblout')
CMscan_df = pd.read_csv(CMscan_output, sep='\t', index_col=None, low_memory=False, usecols=[0,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17], header=None, names=["gene", "accession","model", "mdl_from", "mdl_to", "seq_from", "seq_to", "strand", "trunc", "pass", "gc", "bias", "score", "E-value", "Inc", "Length"])
CMscan_df = CMscan_df[1:]
#print(CMscan_df.head(10))

#Remove 5S model rows
CMscan_df2 = CMscan_df[CMscan_df['gene'] != "5S_rRNA"]
#print(CMscan_df2.head(20))
#Find sequences on the minus strand
minus_strand = CMscan_df2[CMscan_df2['strand'] != "+"]
print("I found sequences on the minus strand ", minus_strand)
#NEEDS add a step to flip any sequences found on the minus strain

#Find sequences that have truncated models
truncated = CMscan_df2[CMscan_df2['trunc'] == "5'&3'"]
print("I found truncated models suggesting the presence of an intron ", truncated)
#print(truncated['accession'])
#Possible action needed here

#Find sequences that do not pass CMscan tests
fail_test = CMscan_df2[CMscan_df2['Inc'] == "?"]
print("Sequences that have a ? are ", fail_test)
#Possible action needed here to exclude some sequences

#show rows containing SSU or LSU
SSU_RNA_df = CMscan_df2[CMscan_df2['gene'] == "SSU_rRNA_eukarya"]
#print("these sequences have SSU ", SSU_RNA_df)

Five_RNA_df = CMscan_df2[CMscan_df2['gene'] == "5_8S_rRNA"]
#print("these sequences have 5.8S ", Five_RNA_df)

#Sequences with extra sequence on the 5' the extends beyond position 1 of the SSU model
SSU_not5_end = SSU_RNA_df[SSU_RNA_df['seq_from'] != 1]
SSU_model_start = SSU_not5_end[SSU_not5_end['mdl_from'] == 1]
#print(SSU_model_start)
#Sequnces with extra data on the 3' end beyond the end of LSU
LSU_RNA_df = CMscan_df2[CMscan_df2['gene'] == "LSU_rRNA_eukarya"]
LSUextra=LSU_RNA_df.loc[(LSU_RNA_df['seq_to'] != LSU_RNA_df['Length']) & (LSU_RNA_df['mdl_to'] == 3401) & (LSU_RNA_df['mdl_from'] == 1)]

SSU_LSUexact=LSU_RNA_df.loc[(LSU_RNA_df['seq_to'] == LSU_RNA_df['Length']) & (LSU_RNA_df['mdl_to'] == 3401) & (LSU_RNA_df['mdl_from'] == 1)]
SSU_LSUpartial=LSU_RNA_df.loc[(LSU_RNA_df['seq_to'] == LSU_RNA_df['Length']) & (LSU_RNA_df['mdl_to'] < 3401) & (LSU_RNA_df['mdl_from'] > 1)]
#print(LSUextra)

#Trim the sequences with extra sequence flanking the SSU and LSU then rewrite the editted fasta to a new file
from Bio import SeqIO
sequences = []  
for seq_record in SeqIO.parse("ribo-out/ribo-out.ribodbmaker.final.fa", "fasta"): 
    s = seq_record
    if seq_record.id not in Five_RNA_df['accession'].tolist(): 
        #sequences.remove(seq_record)
        print("No 5.8S rRNA was found in ", seq_record.id)
        #missingRNA.append(seq_record)
    else:
        if seq_record.id in LSUextra['accession'].tolist():           
            start_df = CMscan_df2[(CMscan_df2['accession'] == seq_record.id) & (CMscan_df2['gene'] == 'SSU_rRNA_eukarya')]
            start = start_df['seq_from'].iloc[0]
            to_df = CMscan_df2[(CMscan_df2['accession'] == seq_record.id) & (CMscan_df2['gene'] == 'LSU_rRNA_eukarya')]
            to = to_df['seq_to'].iloc[0]            
            s.seq = s.seq[start-1:to]
            #seq_record.description = seq_record.description + " small subunit ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2 and large subunit ribosomal RNA, complete sequence"
            seq_record.description = seq_record.description + " COMPLETE SSU-ITS1-5.8S-ITS2-LSU"
            #print(seq_record.id,seq_record.description)
        else:
            if seq_record.id in SSU_RNA_df['accession'].tolist():  
                seq_record.description = seq_record.description + " SSU"
            if seq_record.id in Five_RNA_df['accession'].tolist():  
                seq_record.description = seq_record.description + " 5.8S"
            if seq_record.id in LSU_RNA_df['accession'].tolist():  
                seq_record.description = seq_record.description + " LSU"
        #else:
            #seq_record.description = seq_record.description + " ITS region"
        sequences.append(s)
SeqIO.write(sequences, outputfile, "fasta")  
#SeqIO.write(missingRNA, "missing5_8.fsa", "fasta")

#Print seq_record.id and sequence length from output file
from Bio import SeqIO
for record in SeqIO.parse(outputfile, "fasta"):  
    tmp = LSUextra[LSUextra['accession'] == record.id]    
    if not tmp.empty:
        length = len(record)
        print("Sequences that were trimmed are: ", record.id, length)
        
#Add the new definition lines to the outputfile by appending text to seq_record.description
        

#Count the number of sequences in the final output file
fh = open(outputfile)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
fh.close()
print("The number of sequences in the final output file is: ", n)


print ("processITS script is done and output is saved in:", outputfile)