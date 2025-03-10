#!\usr\bin\env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 07:59:14 2024

@author: allen
"""
import numpy as np

fileName = '1aki.pdb'

ignoreCards = ['CRYST1', \
               'SEQRES', \
               'SOURCE', \
               'MASTER', \
               'ORIGX3', \
               'AUTHOR', \
               'JRNL', \
               'SCALE1', \
               'TER', \
               'END', \
               'SSBOND', \
               'HEADER', \
               'COMPND', \
               'SCALE2', \
               'KEYWDS', \
               'HETATM', \
               'ORIGX1', \
               'EXPDTA', \
               'FORMUL', \
               'REMARK', \
               'HELIX', \
               'SHEET', \
               'REVDAT', \
               'DBREF', \
               'CONECT', \
               'ORIGX2', \
               'TITLE', \
               'SCALE3' \
               ]

coords = np.loadtxt( fileName, \
                     comments = ignoreCards, \
                     usecols=(6,7,8) 
                     )


        
