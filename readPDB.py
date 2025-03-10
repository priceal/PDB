#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 07:59:14 2024

@author: allen
"""


fileName = 'pdb/1aki.pdb'
select = 'ca'
chain = 'A'

with open(fileName,'r') as file:
    lines = file.readlines()
    
match select:
    case 'all':
        
        coordinates = []  
        for line in lines:
            lineSplit = line.split()
            match lineSplit[0]:
                case 'ATOM':
                    chainId = line[21:22]
                    if chainId.strip()==chain:
                        coordinates.append([line[30:38],line[38:46],line[46:54]])
                case _:
                    pass
        result = np.array(coordinates,dtype=float)
        
    case 'ca':
        
        coordinates = []  
        for line in lines:
            lineSplit = line.split()
            match lineSplit[0]:
                case 'ATOM':
                    atomName = line[12:16]
                    chainId = line[21:22]
                    if atomName.strip() == 'CA' and chainId.strip()==chain:
                        coordinates.append([line[30:38],line[38:46],line[46:54]])
                case _:
                    pass
        result = np.array(coordinates,dtype=float)
        
    case 'sequence':
        
        sequence = []
        for line in lines:
            lineSplit = line.split()
            match lineSplit[0]:
                case 'ATOM':
                    atomName = line[12:16]
                    resName = line[17:20]
                    chainId = line[21:22]
                    if atomName.strip() == 'CA' and chainId.strip()==chain:
                        sequence.append(resName)
                case _:
                    pass
        result = sequence


    case 'test':
       
        sequence = []
        for line in lines:
            lineSplit = line.split()
            match lineSplit[0]:
                case 'ATOM':
                    atomName = line[12:16]
                    resName = line[17:20]
                    chainId = line[21:22]
                    seqNumber = int(line[22:26])
                    xyz=[float(x) for x in (line[30:38],line[38:46],line[46:54])]
                    coordinates.append(xyz)
                    if atomName.strip() == 'CA':
                        sequence.append(resName)
                case _:
                     pass
        result = 'test'