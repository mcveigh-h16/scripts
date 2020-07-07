# -*- coding: utf-8 -*-
"""
Created on Fri May 22 13:08:02 2020

@author: mcveigh
"""

import pandas as pd
import os
import sys

inputfile = sys.argv[1]
outputfile = sys.argv[2]

accession = []

clusteruc_df = pd.read_csv(inputfile, sep='\t', index_col=None, low_memory=False, header=None, usecols=[0,2,8], names=["type","cluster","accession"])
clusteruc_df['accession'] = clusteruc_df['accession'].astype(str)
clusteruc_df['cluster'] = pd.to_numeric(clusteruc_df['cluster'])

#print(clusteruc_df)
accession = clusteruc_df['accession']
#accession.to_string()
#print(accession)

with open('acclist', 'w') as filehandle:
    for listitem in accession:
        filehandle.write('%s\n' % listitem)
os.system("/netopt/ncbi_tools64/bin/srcchk -i acclist -f taxname -o acclist.taxdata")

taxdata_file_name = (r'acclist.taxdata')    
taxdata_df = pd.read_csv(taxdata_file_name, sep='\t', index_col=None, header=None, low_memory=False, skiprows=1, usecols=[0,1], names=["accession","taxname"])
#header = taxdata_df.iloc[0]
#taxdata_df.rename(columns = header)

#taxdata_df.to_csv('acclist.taxmap', sep='\t', index=False, header=False, columns=['accession','taxid'])
#print(taxdata_df)
merged_df = pd.merge(left=clusteruc_df, right=taxdata_df, left_on='accession', right_on='accession')
#print(merged_df)
merged_df.to_csv(outputfile, sep='\t', index=False, header=True)

#count = merged_df.groupby('taxname').count()
#merged_df['cluster'] = pd.to_numeric(merged_df['cluster'])
singletons = merged_df[merged_df['cluster'] == 1]
#print(singletons)

multi_df = merged_df[merged_df['cluster'] > 1]

single_multi = pd.merge(singletons, multi_df, on=['taxname'], how='inner')
print(single_multi)
single_multi.to_csv(r'single_multi_intersec' , sep='\t', index=False, header=True)
        