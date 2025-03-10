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

if not os.path.exists(dataDirectory): exit(1)

xyz = []   
for code in pdbCodes:
    path=os.path.join(dataDirectory,code+'.csv')
#    xyz.append(np.load(path))
    df = pd.read_csv(path)
    
    xyz = df[['x','y','z']].to_numpy()
    pt.plotCoords(xyz,ltype='.')