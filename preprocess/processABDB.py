#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:44:14 2025

downloads PDBs and creates (if needed) directory to store

@author: allen
"""

summaryFile = 'final_ab-ag_summary_short.csv'
dataDirectory = 'ABDB_csv'

######################################################################
import os
import numpy as np
import pandas as pd
import pdbTools as pt
######################################################################

summary = pd.read_csv(summaryFile)
pdbCodes = [ cd.lower() for cd in summary['PDB_Code'][:100] ]

if not os.path.exists(dataDirectory):
    os.mkdir(dataDirectory)
    
for code in pdbCodes:
    data = pt.load(code)
#    xyz = pt.getData(data, select='ca')
    df = pt.makeDataFrame( data )
    path=os.path.join(dataDirectory,code+'.csv')
#    np.save(path,xyz)
    df.to_csv(path)
   