#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:44:14 2025

@author: allen
"""

pdbCodes = ['1qc9','6ba1']
dataDirectory = 'testdir'

######################################################################
import os
import numpy as np
import pdbTools as pdb
######################################################################
data = {}
xyz = {}

os.mkdir(dataDirectory)
for code in pdbCodes:
    data[code] = pdb.fetch(code)
    xyz[code] = pdb.extractCoords( data[code], ca=True)
    path=os.path.join(dataDirectory,code)
    np.save(path,xyz[code])


