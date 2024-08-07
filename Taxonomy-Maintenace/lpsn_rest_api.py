# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 09:58:41 2024

@author: mcveigh

LPSN rest API download all new entries by date, compare to NCBI taxonomy with 
gettax. Merge output an highlight names missing from NCBI taxonomy
"""

import pandas as pd
import lpsn
import requests
import re
import json
from pandas import json_normalize
import os
import sys
import subprocess
import xlsxwriter

pub_df = pd.DataFrame(columns=['tid', 'monomial', 'species_epithet', 'subspecies_epithet', 'full_name', 'authority', 'category',
                              'proposed_as', 'validly_published', 'is_legitimate', 'nomenclatural_status', 'is_spelling_corrected',
                              'nomenclatural_type_id', 'type_strain_names', 'basonym_id', 'publication_text', 'publication_kind',
                              'ijsem_list_text', 'ijsem_list_kind', 'emendations', 'molecules', 'lpsn_correct_name_id',
                              'lpsn_taxonomic_status', 'lpsn_address'])

pd.set_option('display.max_columns', None)

client = lpsn.LpsnClient('mcveigh@ncbi.nlm.nih.gov', 'Ysgrtb`1')
result = []
# the prepare method fetches all LPSN-IDs matching your query
# and returns the number of IDs found
#count = client.search(taxon_name='Sulfolobus', correct_name='yes')
count = client.search(date_option='after', date='2024-06-01', correct_name='yes')
print(count, 'entries found.')

# the retrieve method lets you iterate over all results
# and returns the full entry as dict
# Entries can be further filtered using a list of keys (e.g. ['keywords'])
for entry in client.retrieve():
    #print(entry)
    result.append(entry)
#print(result)

df = json_normalize(result)
lpsn_df = df.drop(columns=['id', 'molecules', 'species_epithet', 'lpsn_correct_name_id', 'nomenclatural_type_id', 'basonym_id', 'subspecies_epithet'])
lpsn_df.rename(columns={'full_name' : 'lpsn_name'}, inplace=True)
#print(lpsn_df.head)

name_list = lpsn_df['lpsn_name'].values.tolist()
#print(name_list)

gettax_df = pd.DataFrame(columns=['lpsn_name', 'ncbi_name'])

for lpsn_name in name_list:
    ncbi_name = []
    #print('name is', name)
    args = 'gettax -1 ' + lpsn_name
    p = subprocess.Popen(args, shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    out, err = p.communicate()
    gettax_response = str(out)
  
    if gettax_response == "":
        print (lpsn_name, 'is not found in NCBI taxonomy')
        ncbi_name = "not found"
    else:
        sub1 = ":"
        sub2 = "|"
        idx1 = gettax_response.index(sub1)
        idx2 = gettax_response.index(sub2)
        ncbi_name = gettax_response[idx1 + len(sub1) + 1: idx2]
    #print('names found', lpsn_name, ncbi_name)
    #print ('error', err)
    
    row_data = [lpsn_name, ncbi_name]
    length = len(gettax_df)
    gettax_df.loc[length] = row_data

#print(gettax_df.head)

#merge dataframes
combine_df=pd.merge(left=lpsn_df, right=gettax_df, left_on='lpsn_name', right_on='lpsn_name', how = 'outer')
combine_df = combine_df[['lpsn_name', 'ncbi_name', 'category', 'monomial', 'authority', 'emendations', 'proposed_as', 'lpsn_address', 'is_legitimate', 
'ijsem_list_kind', 'ijsem_list_text', 'publication_kind', 'publication_text', 'type_strain_names',  'validly_published', 'nomenclatural_status',
'is_spelling_corrected', 'lpsn_taxonomic_status']]
#print(combine_df.head)

def highlight_rows(row):
    lspnvalue = row.loc['lpsn_name']
    ncbivalue = row.loc['ncbi_name']
    if lspnvalue != ncbivalue:
        color = '#FFB3BA' # Red
    elif lspnvalue == ncbivalue:
        color = '#BAFFC9' # Green
    return ['background-color: {}'.format(color) for r in row]

new_df = combine_df.style.apply(highlight_rows, axis=1, subset=['lpsn_name', 'ncbi_name'])

#lpsn_df.to_excel("lpsn_output.xlsx") 
#combine_df.to_excel("lpsncombined_output.xlsx") 
new_df.to_excel('lpsncombinestyle.xlsx', engine='xlsxwriter', index = False, na_rep = '')     
    