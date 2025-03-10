#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:44:14 2025

downloads PDBs and creates (if needed) directory to store

@author: allen
"""

pdbCodes = ['1qc9','6ba1','2lyz']
dataDirectory = 'testdir'

######################################################################
import os
import numpy as np
import pdbTools as pdb
######################################################################
if not os.path.exists(dataDirectory):
    os.mkdir(dataDirectory)
    
for code in pdbCodes:
    data = pdb.fetch(code)
    xyz = pdb.getCoordinates(data, select='ca')
    seq = pdb.getSequence(data)
    saveData = np.concatenate((seq[:,np.newaxis],xyz), axis=1)
    path=os.path.join(dataDirectory,code)
    np.save(path,saveData)
   