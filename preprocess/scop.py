#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 08:44:59 2025

@author: allen
"""

scopFile = 'data/scop-cla-latest.csv'
sampleSize = 1



##########################################################################
import pandas as pd

scop=pd.read_csv(scopFile)

scopSample = scop.sample( sampleSize )
scopSet = set(scopSample['FA-PDBID'])
print(scopSet)
print(len(scopSet))

