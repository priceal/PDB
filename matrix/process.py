#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:44:14 2025

downloads PDBs and creates (if needed) directory to store

@author: allen
"""

dataDirectory = '/home/allen/projects/PDB/preprocess/scop10_csv'

######################################################################
import os
import sys
import numpy as np
import pandas as pd
import pdbTools as pt
######################################################################
if not os.path.exists(dataDirectory):
    print(dataDirectory+' does not exist')
    sys.exit()

fileNames = os.listdir(dataDirectory)
paths = [ os.path.join(dataDirectory, names) for names in fileNames ]

aaContacts = np.zeros((20,20))   
for i, path in enumerate(paths):
    print(str(i)+'. '+path)
    df = pd.read_csv(path)
    dfca = df[ df['name']=='CA' ]
    xyz = dfca[ ['x','y','z'] ].to_numpy()
    seqRaw = dfca[ 'res' ].to_list()
    seq = [ pt.defaultCode[ res ] for res in seqRaw ]
    aaContacts += pt.aaContactMatrix(xyz, seq)

#now fix the central diagonal
aaContacts = aaContacts - 0.5*np.diag(aaContacts.diagonal())

pt.heatMap(aaContacts,tickLabels=pt.defaultCode.keys())  

