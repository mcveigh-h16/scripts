# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 14:41:14 2019

@author: mcveigh
"""

import pandas as pd
import numpy as np

file_name_string = (r'C:\Users\mcveigh\Documents\PythonPC\All_ITS_tax_report.xlsx')
taxreport_df = pd.read_excel(file_name_string, index_col=None, na_values=['-'])

taxreport_df.groupby('kingdom').count()