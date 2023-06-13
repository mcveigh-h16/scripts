# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 08:59:25 2023

@author: mcveigh

revisision of processITS1.98.py to better handle reverse complemented sequences
any sequences with 5.8S rRNA on the minus strand are reverse complemented. Rerun cmscan on these
and merge the results. Requires strong match to complete 5.8S rRNA for inclusion. 

"""
import pandas as pd
import Bio
import os
import sys
from datetime import datetime
import functools

startTime = datetime.now()
print("Start time is ", startTime) 

inputfile = sys.argv[1]
outputfile = sys.argv[2]
#inputfile = (r'test2.gbk')
#outputfile = (r'test.out')
"""
Read in the reject list and find the accession
"""
reject_file_name = (r'ITS_reject_seqs3.txt') 
rejectlist_df = pd.read_csv(reject_file_name, sep='\t', index_col=None, low_memory=False, header=None, names=["accession", "type", "reason"])
rejectlist = rejectlist_df['accession']
reject_list = set(rejectlist_df['accession'].tolist())

"""   
Parse the GenBank file and remove any sequences found on the reject list
"""
from Bio import SeqIO
sequences = [] 
found = []
sequencelength = []

rejectTime = datetime.now()
print("Reject time is ", rejectTime) 

for seq_record in SeqIO.parse(inputfile, "genbank"):    
    str_id = seq_record.id     
    if seq_record.name not in reject_list:                       
        seq_record.description = seq_record.annotations["organism"]
        sequences.append(seq_record)
        seqlength = '%s %i\n' %  (seq_record.id, len(seq_record))
        sequencelength.append(seqlength)
    else:
        found.append(seq_record)
        print("I found this accession on the reject list and wrote the sequence to found.fsa: ", seq_record.id)  
SeqIO.write(sequences, "stripped.fsa", "fasta")    
seqlen_str = functools.reduce(lambda a,b : a + b, sequencelength) 
f = open('my.seqlen', 'w')
f.write(seqlen_str)
f.close()
  
"""
Run short version of ribodbmaker.pl checks for specified species, internal Ns and vector contamination   
"""
os.system("ribodbmaker -f --skipfribo1 --skipfribo2 --skipfmspan --skipingrup --skipclustr --skiplistms --skipmstbl --skipfblast --taxin /panfs/pan1/dnaorg/rrna/git-ncbi-rrna-project/taxonomy-files/ncbi-taxonomy-tree.ribodbmaker.txt stripped.fsa ribo-out")
#os.system("ribodbmaker -f --skipfribo1 --skipfribo2 --skipfmspan --skipingrup --skipclustr --skiplistms --skipmstbl --skipfblast --skipftaxid  stripped.fsa ribo-out")
#os.system("ribodbmaker -f --skipfribo1 --skipfribo2 --skipfmspan --skipingrup --skipclustr --skiplistms --skipmstbl --skipfblast stripped.fsa ribo-out")
riboTime = datetime.now()
print("ribodbmaker time is ", riboTime) 

"""
Run CMscan using the output of ribodbmaker.pl
"""
print('cmscan1')
os.system("cmscan --cpu 16 --mid -T 20 --verbose --tblout tblout.df.txt rrna.cm ribo-out/ribo-out.ribodbmaker.full.fa > /dev/null")
print('cmscan2')
os.system("cmscan --cpu 16 --mid -T 20 --verbose --anytrunc --tblout tblout.at.txt rrna.cm ribo-out/ribo-out.ribodbmaker.full.fa > /dev/null")
cmscanTime = datetime.now()
print("CMscan time is ", cmscanTime) 
os.system("cat tblout.df.txt tblout.at.txt > tblout.both.txt")
os.system("perl cmsearch_tblout_deoverlap/cmsearch-deoverlap.pl --maxkeep -s --cmscan tblout.both.txt")
os.system("head -n2 tblout.both.txt > final.tblout")
os.system("cat tblout.both.txt.deoverlapped >> final.tblout")

"""
Add seq length to cmscan output
"""
os.system("perl tblout-add.pl -t final.tblout 18 my.seqlen 3 > cmscan_final.tblout")

"""
Parse the results of CMscan. If any 5.8S rRNA are found on minus strand these
are reverse complemented and rerun through CMscan
"""
CMscan_output = (r'cmscan_final.tblout')
CMscan_df = pd.read_csv(CMscan_output,
                        sep='\t',
                        index_col=None,
                        low_memory=False,
                        usecols=[0,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17],
                        header=None,
                        names=["gene", "accession","model", "mdl_from",
                               "mdl_to", "seq_from", "seq_to", "strand",
                               "trunc", "pass", "gc", "bias", "score",
                               "E-value", "Inc", "Length"])


FiveCompleteMinus = CMscan_df.loc[(CMscan_df['gene'] == "5_8S_rRNA") & (CMscan_df['trunc'] == "no") & (CMscan_df['score'] > 50) & (CMscan_df['strand'] != "+")]
#print(FiveCompleteMinus)

reverse_seq = []
reverse_acc = []
forward_seq = []

for seq_record in SeqIO.parse("ribo-out/ribo-out.ribodbmaker.full.fa", "fasta"): 
    s = seq_record
    if seq_record.id in FiveCompleteMinus['accession'].tolist():
        print("I reverse complemented ", seq_record.id)
        s.seq = s.seq.reverse_complement()
        reverse_seq.append(s)
        reverse_acc.append(s.id)
    else:
        forward_seq.append(s)
        
SeqIO.write(forward_seq, "forward_seqs", "fasta")
SeqIO.write(reverse_seq, "reverse_seqs", "fasta")
os.system("cat forward_seqs reverse_seqs > plus_strand_seqs")

def repeat_cmscan():
    """
    repeat cmscan for only sequences that have been reverse complemented and
    combines output original. 
    """  

    CMscan_df.drop(CMscan_df[CMscan_df.accession.isin(reverse_acc)].index.tolist())
    Plus_df = CMscan_df[~CMscan_df.accession.isin(reverse_acc)]
    

    """
    Run CMscan on the reverse complemented sequences
    """
    print('cmscan3')
    os.system("cmscan --cpu 16 --mid -T 20 --verbose --tblout repeat.df.txt rrna.cm reverse_seqs > /dev/null")
    print('cmscan4')
    os.system("cmscan --cpu 16 --mid -T 20 --verbose --anytrunc --tblout repeat.at.txt rrna.cm reverse_seqs > /dev/null")
    repeatcmscanTime = datetime.now()
    print("CMscan time is ", repeatcmscanTime) 
    os.system("cat repeat.df.txt repeat.at.txt > repeat.both.txt")
    os.system("perl cmsearch_tblout_deoverlap/cmsearch-deoverlap.pl --maxkeep -s --cmscan repeat.both.txt")
    os.system("head -n2 repeat.both.txt > repeat.tblout")
    os.system("cat repeat.both.txt.deoverlapped >> repeat.tblout")
    os.system("perl tblout-add.pl -t repeat.tblout 18 my.seqlen 3 > repeat_cmscan_final.tblout")

    """
    Collect reverse complemented CMscan data and merging this with the data for
    sequences on the plus strand. All seqs should now be plus strand 5.8S rRNA
    """
    repeat_CMscan_output = (r'repeat_cmscan_final.tblout')
    repeat_CMscan_df = pd.read_csv(repeat_CMscan_output,
                        sep='\t',
                        index_col=None,
                        low_memory=False,
                        usecols=[0,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17],
                        header=None,
                        names=["gene", "accession","model", "mdl_from",
                               "mdl_to", "seq_from", "seq_to", "strand",
                               "trunc", "pass", "gc", "bias", "score",
                               "E-value", "Inc", "Length"])

    combined = [Plus_df, repeat_CMscan_df]
    combined_df = pd.concat(combined)
    print(combined_df)
    print('i finished the function')
    return combined_df
        
if len(reverse_acc) != 0:
    #print('Minus strand sequences are found')
    print('reversed accesssions are: ', reverse_acc)
    repeat_cmscan()
    combined_df = repeat_cmscan()
else:
    combined_df = CMscan_df
    print('creating combinedf')

#CREATE FUNCTION FOR MINUS STRAND was here

"""
Remove 5S model rows from data frame
Find and report sequences that have truncated models
Find sequences that do not pass CMscan tests
"""

combined_df = combined_df[combined_df['gene'] != "5S_rRNA"]
print(combined_df)

truncated=combined_df.loc[(combined_df['gene'] != "5_8S_rRNA") & (combined_df['trunc'] == "5'&3'")]
print("I found truncated models suggesting the presence of an intron\n ", truncated)

fail_test = combined_df[combined_df['Inc'] == "?"]
print("Sequences that have a ? are\n ", fail_test)

"""
Build dataframes for further use parsing sequences
"""
Fivecomplete = combined_df.loc[(combined_df['gene'] == "5_8S_rRNA") & (combined_df['trunc'] == "no") & (combined_df['score'] > 50)]
SSU_RNA_df = combined_df.loc[(combined_df['gene'] == "SSU_rRNA_eukarya")]
LSU_RNA_df = combined_df.loc[(combined_df['gene'] == "LSU_rRNA_eukarya")]  
LSUpartial = combined_df.loc[(combined_df['gene'] == "LSU_rRNA_eukarya") & (combined_df['trunc'] != "no")]
LSUcomplete = combined_df.loc[(combined_df['gene'] == "LSU_rRNA_eukarya") & (combined_df['trunc'] == "no")]
SSUpartial = combined_df.loc[(combined_df['gene'] == "SSU_rRNA_eukarya") & (combined_df['trunc'] != "no")]
SSUcomplete = combined_df.loc[(combined_df['gene'] == "SSU_rRNA_eukarya") & (combined_df['trunc'] == "no")]                         
LSUendmissing = LSUpartial[(LSUpartial['strand'] == "+") & (LSUpartial['seq_to'] != LSUpartial['Length'])]
SSUendmissing = SSUpartial[(SSUpartial['strand'] == "+") & (SSUpartial['seq_from'] > 1)]
SSUendfound = SSUpartial[(SSUpartial['strand'] == "+") & (SSUpartial['mdl_from'] == 1)]
LSUendfound = LSUpartial[(LSUpartial['strand'] == "+") & (LSUpartial['mdl_to'] == 3401)]
SSUminus = combined_df.loc[(combined_df['gene'] == "SSU_rRNA_eukarya") & (combined_df['strand'] != "+")]
LSUminus = combined_df.loc[(combined_df['gene'] == "LSU_rRNA_eukarya") & (combined_df['strand'] != "+")]
SSUextra = SSUcomplete.loc[(SSUcomplete['seq_from'] != 1) & SSUcomplete['mdl_from'] == 1]
LSUextra=LSUcomplete.loc[(LSUcomplete['seq_to'] != LSUcomplete['Length']) & (LSUcomplete['mdl_to'] == 3401) & (LSUcomplete['mdl_from'] == 1)]
print("sequences with extra data on the 5' end \n", SSUextra)
print("sequences with extra data on the 3' end \n", LSUextra)


seqdata = []  
misassembled = []
removeacc = []
#rna_not_found = []
rna_found = []

"""
Parse all sequences, look for misassemblies, trim extra sequences on ends and build definition lines
"""


for seq_record in SeqIO.parse("plus_strand_seqs", "fasta"): 
    s = seq_record
    SSUstart = []
    SSUend = []
    fivestart = []
    fiveend = []
    LSUstart = []
    LSUend = []
    
    if seq_record.id not in Fivecomplete['accession'].tolist(): 
        #sequences.remove(seq_record)
        print("No complete 5.8S rRNA was found in", seq_record.id)

    """            
    Check for Mixed Strand, Noncontiguous and Misassembled Sequences 
    """
    if seq_record.id in Fivecomplete['accession'].tolist():
        if seq_record.id in SSUminus['accession'].tolist():
            print(seq_record.id, "Mixed strand sequence 1")
            misassembled.append(s)
            removeacc.append(seq_record.id)
        if seq_record.id in LSUminus['accession'].tolist():
            print(seq_record.id, "Mixed strand sequence 2")
            misassembled.append(s)
            removeacc.append(seq_record.id)
        if seq_record.id in SSU_RNA_df['accession'].tolist():
            start_df = SSU_RNA_df[(SSU_RNA_df['accession'] == seq_record.id) & (SSU_RNA_df['gene'] == 'SSU_rRNA_eukarya')]
            SSUstart = int(start_df['seq_from'].iloc[0])
            SSUend = int(start_df['seq_to'].iloc[0])
            #print(seq_record.id, " SSU ", SSUstart, SSUend)
        if seq_record.id in Fivecomplete['accession'].tolist():
            startfive_df = Fivecomplete[(Fivecomplete['accession'] == seq_record.id) & (Fivecomplete['gene'] == '5_8S_rRNA')]
            fivestart = int(startfive_df['seq_from'].iloc[0])
            fiveend = int(startfive_df['seq_to'].iloc[0])
            #print(seq_record.id, " FIVE ", fivestart, fiveend)
        if seq_record.id in LSU_RNA_df['accession'].tolist(): 
            startLSU_df = LSU_RNA_df[(LSU_RNA_df['accession'] == seq_record.id) & (LSU_RNA_df['gene'] == 'LSU_rRNA_eukarya')]
            LSUstart = int(startLSU_df['seq_from'].iloc[0])
            LSUend = int(startLSU_df['seq_to'].iloc[0])
            #print(seq_record.id, " LSU ", LSUstart, LSUend)
            
        if seq_record.id in SSU_RNA_df['accession'].tolist():
            if seq_record.id in Fivecomplete['accession'].tolist():
                if SSUend > fivestart:
                    print(seq_record.id, " ssu spans ", SSUstart, SSUend)
                    print(seq_record.id, " five spans ", fivestart, fiveend)
                    print(seq_record.id, "Misassembled sequence 1")
                    misassembled.append(s)
                    removeacc.append(seq_record.id)
        if seq_record.id in Fivecomplete['accession'].tolist():
            if seq_record.id in LSU_RNA_df['accession'].tolist():
                if fiveend > LSUstart:
                    print(seq_record.id, "Misassembled sequence 2")
                    print(seq_record.id, " lsu spans ", LSUstart, LSUend)
                    misassembled.append(s)
                    removeacc.append(seq_record.id)   
    SeqIO.write(misassembled, "misassembled_seqs", "fasta") 
    if seq_record.id not in removeacc:
        rna_found.append(s)
        SeqIO.write(rna_found, "rna_found.seqs", "fasta")   
        
    """
    Trim extra sequence upstream of SSU or downstream of LSU
    """        
    if seq_record.id in Fivecomplete['accession'].tolist():
        if seq_record.id in LSUextra['accession'].tolist():           
            to_df = combined_df[(combined_df['accession'] == seq_record.id) & (combined_df['gene'] == 'LSU_rRNA_eukarya')]
            to = to_df['seq_to'].iloc[0] 
            s.seq = s.seq[0:int(to)]

        if seq_record.id in SSUextra['accession'].tolist():
            start_df = SSUextra[(SSUextra['accession'] == seq_record.id) & (SSUextra['gene'] == 'SSU_rRNA_eukarya')]
            start = start_df['seq_from'].iloc[0]
            to = start_df['Length'].iloc[0]
            s.seq = s.seq[start-1:to]
  
    if seq_record.id in LSUendmissing['accession'].tolist(): 
        if seq_record.id not in LSUendfound['accession'].tolist(): 
            to_df = combined_df[(combined_df['accession'] == seq_record.id) & (combined_df['gene'] == 'LSU_rRNA_eukarya')]
            to = to_df['seq_to'].iloc[0]
            s.seq = s.seq[0:int(to)]
                
    if seq_record.id in SSUendmissing['accession'].tolist():  
        if seq_record.id not in SSUendfound['accession'].tolist(): 
            if seq_record.id not in SSUextra['accession'].tolist():
                start_df = combined_df[(combined_df['accession'] == seq_record.id) & (combined_df['gene'] == 'SSU_rRNA_eukarya')]
                start = start_df['seq_from'].iloc[0]
                to = start_df['Length'].iloc[0]
                s.seq = s.seq[start-1:to]  
        
    """
    Definition Line generator
    """           
    if seq_record.id not in removeacc:
        if seq_record.id in Fivecomplete['accession'].tolist():
            if seq_record.id in SSUcomplete['accession'].tolist():  
                seq_record.description = seq_record.description + " SSU ITS1 5.8S ITS2"
            elif seq_record.id in SSUendfound['accession'].tolist():   
                seq_record.description = seq_record.description + " SSU ITS1 5.8S ITS2"
            elif seq_record.id in SSUpartial['accession'].tolist(): 
                if seq_record.id not in SSUendfound['accession'].tolist():
                    seq_record.description = seq_record.description + " <SSU ITS1 5.8S ITS2"
            else: 
                seq_record.description = seq_record.description + " <ITS1 5.8S ITS2"       
            if seq_record.id in LSUcomplete['accession'].tolist():  
                seq_record.description = seq_record.description + " LSU"
            elif seq_record.id in LSUpartial['accession'].tolist():  
                seq_record.description = seq_record.description + " LSU>"
            else:
                seq_record.description = seq_record.description + ">"
            seqdata.append(s)
            
SeqIO.write(seqdata, outputfile, "fasta")  
#SeqIO.write(missingRNA, "missing5_8.fsa", "fasta")

"""
Print seq_record.id and sequence length from output file for any sequences that were trimmed
"""
for record in SeqIO.parse(outputfile, "fasta"):  
    tmp = LSUextra[LSUextra['accession'] == record.id]
    tmp2 = SSUextra[SSUextra['accession'] == record.id]
    tmp3 = LSUendmissing[LSUendmissing['accession'] == record.id]
    tmp4 = SSUendmissing[SSUendmissing['accession'] == record.id]
    tmpframe = [tmp, tmp2, tmp3, tmp4]
    extra = pd.concat(tmpframe)
    if not extra.empty:
        length = len(record)
        print("Sequences that were trimmed are: ", record.id, length)
"""        
Count the number of sequences in the final output file
"""
fh = open(outputfile)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
fh.close()
print("The number of sequences in the final output file is: ", n)


print ("processITS script is done and output is saved in:", outputfile)