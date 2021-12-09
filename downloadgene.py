# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 11:42:43 2020

@author: mcveigh
"""

from __future__ import print_function
import time
import ncbi.datasets
from ncbi.datasets.rest import ApiException
from pprint import pprint

# Enter a context with an instance of the API client
with ncbi.datasets.ApiClient() as api_client:
    # Create an instance of the API class
    api_instance = ncbi.datasets.DownloadApi(api_client)
    gene_ids = [56] # list[int] | NCBI Gene ID
include_sequence_type=['SEQ_TYPE_GENE', 'SEQ_TYPE_RNA']
filename = 'filename_example' # str | Output file name. (optional)

try:
# Retrieve a requested gene dataset and stream back reply by gene ID
    api_response = api_instance.download_gene_package(gene_ids, include_sequence_type=include_sequence_type, filename=filename)
    pprint(api_response)
except ApiException as e:
    print("Exception when calling DownloadApi->download_gene_package: %s\n" % e)