#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 08:42:07 2025

@author: allen
"""
import os
import numpy as np

code = '2e2h'

directory = 'data'

fileNames = [file for file in os.listdir(directory) if code in file ]

record={}
for file in fileNames:
    key=os.path.basename(file).split('.')[-2][-1:]
    path =os.path.join(directory,file)
    record[key]=np.load(path)


