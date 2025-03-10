#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:44:14 2025

downloads PDBs and creates (if needed) directory to store

largest side chain Trp: 10 heavy atoms not including backbone (4)
so 14 in all. but could be c-terminal, then 15

@author: allen
"""

pdbCode = '1n8z'
cutoff=4.0

######################################################################
import os
import numpy as np
import pandas as pd
import pdbTools as pt
######################################################################


data = pt.load(pdbCode)
xyz = pt.getData(data, chain='A',select='all')
ca = pt.getData(data,chain='A',select='ca')
pt.plotCoords(xyz,ltype='.')
pt.plotCoords(ca)

df=pt.makeDataFrame(data)

dfA = df[ df['chain'] == 'A' ]
dfAg=dfA.groupby(by='seqnum')

proteinMatrix = np.ones((len(dfAg),15,3))*np.nan
for i, (seqnum, residue) in enumerate(dfAg):
    print( i, seqnum, len(residue) )
#    print( residue[['x','y','z']] )
    proteinMatrix[ i, :len(residue), :] = residue[['x','y','z']].to_numpy()
#    print(proteinMatrix[ i ])

Dv = proteinMatrix[ :         , np.newaxis, :         , np.newaxis, : ] - \
     proteinMatrix[ np.newaxis, :         , np.newaxis, :         , : ]

D  = np.sqrt( (Dv*Dv).sum( axis = 4 ) )

print(D.shape)

Dmin = np.nanmin( D, axis=(2,3))

cutoff=4.0
seq=pt.getData(data,chain='A',select='seq')
# calculate contact matrix & remove self-contacts
CM = np.less(Dmin,cutoff)
np.fill_diagonal(CM, False)

# calculate sequence matrix
seqMatrix = np.zeros( ( len(seq), 20 ) )
seqMatrix[ range(len(seq)), seq] = 1

# perform similarity transform on contact map to xform to aa based map
aaCM = np.linalg.multi_dot([seqMatrix.T,CM,seqMatrix])




  
