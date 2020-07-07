# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 10:14:59 2020

@author: mcveigh
"""

#Script to analyze applog data for accession numbers that have been searched but not found to determine if they should be released

import pandas as pd
import os
import sys
from datetime import datetime
import numpy as np

startTime = datetime.now()
#print("Start time is ", startTime) 

inputfile = sys.argv[1]
#outputfile = sys.argv[2]

starting_df = pd.read_csv(inputfile, sep='\t', index_col=None, low_memory=False, 
                        usecols=[0,4,5], header=None,
                        names=["accession", "insmart","URL"])
#print(starting_df)

accession_df = starting_df['accession']
accession_df.to_string('applog_acclist', index = False, header = False)
os.system("idstat -A applog_acclist -mr -nh -r F -of tbl -i > idstat_report")

idstat_df = pd.read_csv('idstat_report', sep='|', index_col=None, low_memory=False, 
                        usecols=[0,3,4,6,7,8,], header=None,
                        names=["accession", "satellite","status", "withdraw", "suppress", "load"])
#print(idstat_df)
merged_df = pd.merge(left=starting_df, right=idstat_df, left_on='accession', right_on='accession')
merged_df[['loaddate','HUPdate']] = merged_df.load.str.split(", HUP-Date",expand=True)
merged_df = merged_df.drop(columns="load")
merged_df['HUPdate'] = pd.to_datetime(merged_df['HUPdate'])
merged_df['loaddate'] = pd.to_datetime(merged_df['loaddate'])        
merged_df['HUPdaysLeft'] = (merged_df['HUPdate'] - merged_df['loaddate']).dt.days      

#merged_df[['loaddate','HUPdate']] = merged_df.load.apply(lambda x: pd.Series(str(x).split(", HUP-Date"))) 
#merged_df.to_csv(outputfile, sep='\t', index=False, header=True)
#print("The merged data frame is \n", merged_df)

df_all = merged_df.merge(starting_df.drop_duplicates(), on=['accession', 'URL', 'insmart'], how='outer', indicator=False)

#df_all = df_all.drop(df_all[df_all.satellite == "DDBJ_WGS"].index)
#df_all = df_all.drop(df_all[df_all.satellite == "WGS2"].index)
df_all = df_all[~df_all.satellite.str.contains("DDBJ", na=False)]
df_all = df_all[~df_all.satellite.str.contains("WGS", na=False)]
df_all = df_all[~df_all.satellite.str.contains("EMBL", na=False)]
wgs_accession = "(^[A-Za-z]{4}\d{8})"
genome_accession = '^CP(\d{6})'
df_all = df_all[~df_all.accession.str.contains(wgs_accession, na=False)]
df_all = df_all[~df_all.accession.str.contains(genome_accession, na=False)]
#df_all['URL'] = df_all['URL'].str.replace('(?:https?://|www\d{0,3}[.])boldsystem', '')


#urldrop = ['boldsystem', 'barcode', 'seaphage', 'google']
#df_all = df_all[~df_all.URL.str.contains("boldsystem", na=False)]
#df_all = df_all[~df_all.URL.str.contains("barcode", na=False)]
#df_all = df_all[~df_all.URL.str.contains("google", na=False)]
#df_all = df_all[~df_all.URL.str.contains("seaphage", na=False)]

#print(df_all.head)


#df_all.to_string('testoutput', index = False, header = False)
writer = pd.ExcelWriter('AppLog_output.xlsx')
df_all.to_excel(writer,'Sheet1')
writer.save()
#merged_df.to_excel(r'merged.xlsx', index = False)
#df_all.to_excel(r'merged_all.xlsx', index = False)
print('ZeroSearch data is written successfully to Excel file AppLog_output.xlsx')