#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:46:56 2024

@author: allen
"""

import requests
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

defaultCode = { 'ALA': 0, 
                'ARG': 1,
                'ASN': 2,
                'ASP': 3,
                'CYS': 4,
                'GLN': 5,
                'GLU': 6,
                'GLY': 7,
                'HIS': 8,
                'ILE': 9,
                'LEU': 10,
                'LYS': 11,
                'MET': 12,
                'PHE': 13,
                'PRO': 14,
                'SER': 15,
                'THR': 16,
                'TRP': 17,
                'TYR': 18,
                'VAL': 19
            }

### 
propCode = { 'ALA': 0, 
                'ARG': 1,
                'ASN': 2,
                'ASP': 3,
                'CYS': 4,
                'GLN': 5,
                'GLU': 6,
                'GLY': 7,
                'HIS': 8,
                'ILE': 2,
                'LEU': 3,
                'LYS': 11,
                'MET': 4,
                'PHE': 13,
                'PRO': 14,
                'SER': 15,
                'THR': 16,
                'TRP': 17,
                'TYR': 18,
                'VAL': 1
            }







###############################################################################
def fetchSequence( source ):
    '''
    downloads FASTA file from PDB and returns as string

    Args:
        source (str): the 4 character PDB code 

    Returns:
        text (str): the text of the file (or 'error').
        
    '''
    text ='Error'  # define in case return after an error  
    url = 'https://www.rcsb.org/fasta/entry/'+source
    try:
        download = requests.get(url)
    except Exception as exc:
        print(f'exception dowloading \'{url}\'\nException: {exc}')
        return 'Exception'
    else:
        if download.reason=='OK':
            text = download.text
        else:
            print(f'returned status not OK downloading \'{url}\'\nStatus: {download.reason}')
            return 'Bad status'
    return text

#####################################################



















###############################################################################
def load( source, ba='' ):
    '''
    loads PDB file and returns a list of strings, one entry
    per line. source can be either a pdb 4-letter code (in which case the
    file is downloaded) or a path to a pdb file. paths must end with '.pdb'

    Args:
        source (str): the 4 character PDB code or path to file

    Returns:
        lines (list): the text of the file (empty on error).
        
    '''
    lines = [] # define in case return after an error
    # if filename given, read in file and store each line in lines
    if source.endswith('.pdb'):  
        if os.path.exists(source):
            with open(source,'r') as file:
                lines = file.readlines()
        else:
            print('error: file does not exist')
            return lines
    else: 
        url = 'https://files.rcsb.org/download/'+source+'.pdb'+ba
        try:
            download = requests.get(url)
        except Exception as exc:
            print(f'exception while attempting to download \'{url}\' :\n{exc}')
            return lines
        else:
            if download.reason=='OK':
                lines = download.text.split('\n')
            else:
                print('returned status not OK while attempting to download: '+url)
                return lines
    
    return lines

###############################################################################
def getData( source, chain='A', select='all'):
    '''
    extracts coordinates or sequence from source. 
    
    Args:
        source (str or list): path to the PDB file or a list containing
            the lines of the PDB file
        chain (str, optional): which chain to read. Defaults to 'A'.
        select (string, optional): choices are 'all', 'ca', 'seq', 'ca,seq'. 
            Defaults to 'all'.

    Returns:
        numpy array of coordinates, shape = (N,3) or
        list of sequence 3-letter codes, shape = (N) or
        coordinates, sequence: ordered tuple of both

    '''   
    result = [] # define in case return after an error
    # if filename given, read in file and store each line in lines
    if type(source) is str:
        if os.path.exists(source):
            with open(source,'r') as file:
                lines = file.readlines()
        else:
            print('error: file does not exist')
            return result
            
    # else assume the lines of the data are given
    else: lines=source   
    
    # now do the separate cases
    match select:
        case 'all': # store all coordinates in chain
            
            coordinates = []  
            for line in lines:
                if line == '': continue
                match line.split()[0].strip():
                    case 'ATOM':
                        chainId = line[21:22]
                        if chainId.strip()==chain:
                            coordinates.append([line[30:38],line[38:46],line[46:54]])
                    case _:
                        pass
            result = np.array(coordinates,dtype=float)
            
        case 'ca': # only store CA coordinates
            
            coordinates = []  
            for line in lines:
                if line == '': continue
                match line.split()[0].strip():
                    case 'ATOM':
                        atomName = line[12:16]
                        chainId = line[21:22]
                        if atomName.strip() == 'CA' and chainId.strip()==chain:
                            coordinates.append([line[30:38],line[38:46],line[46:54]])
                    case _:
                        pass
            result = np.array(coordinates,dtype=float)
            
        case 'seq': # if sequence is requested, store sequence codes
            
            sequence = []
            for line in lines:
                if line == '': continue
                match line.split()[0].strip():
                    case 'ATOM':
                        atomName = line[12:16]
                        resName = line[17:20]
                        chainId = line[21:22]
                        if atomName.strip() == 'CA' and chainId.strip()==chain:
                            sequence.append( defaultCode[ resName ] )  
                    case _:
                        pass
            result = sequence
   
        case 'ca,seq': # store CA coordinates and sequence
            
            coordinates = []
            sequence = []
            for line in lines:
                if line == '': continue
                match line.split()[0].strip():
                    case 'ATOM':
                        atomName = line[12:16]
                        resName = line[17:20]
                        chainId = line[21:22]
                        if atomName.strip() == 'CA' and chainId.strip()==chain:
                            coordinates.append([line[30:38],line[38:46],line[46:54]])
                            sequence.append(resName)
                    case _:
                        pass
            result = np.array(coordinates,dtype=float), sequence
    
    return result

###############################################################################
def makeDataFrame( source ):
    
    columns = {'record': (1,4), 
               'number': (7,11),
               'name': (13,16),
               'alt': (17,17),
               'res': (18,20),
               'chain': (22,22),
               'seqnum': (23,26),
               'insert': (27,27),
               'x': (31,38),
               'y': (39,46),
               'z': (47,54),
               'occ': (55,60),
               'B': (61,66),
               'segid': (73,76),
               'element': (77,78),
               'charge': (79,80)
                }
    
    dataTypes = {'record': str, 
               'number': int,
               'name': str,
               'alt': str,
               'res': str,
               'chain': str,
               'seqnum': int,
               'insert': str,
               'x': float,
               'y': float,
               'z': float,
               'occ': float,
               'B': float,
               'segid': str,
               'element': str,
               'charge': str
                }
    
    data = dict( zip( columns.keys(), [ [] for i in range(len(columns)) ] ) )
    
    result = [] # define in case return after an error
    
    # if filename given, read in file and store each line in lines
    if type(source) is str:
        if os.path.exists(source):
            with open(source,'r') as file:
                lines = file.readlines()
        else:
            print('error: file does not exist')
            return result
            
    # else assume the lines of the data are given
    else: lines=source   
    
    for line in lines:
        if line == '': continue
        if line.split()[0].strip()=='ATOM':
            for item in columns:
                data[item].append(line[columns[item][0]-1:columns[item][1]].strip())

    return pd.DataFrame( data ).astype( dataTypes )
        
###############################################################################
def plotCoords( coords, ltype='-' ):
    '''
    creates 3D chain trace of coordinates

    Args:
        coords (array (N,3)): the coordinates
        ltype (str, optional): line style. Defaults to '-'.

    Returns:
        None.

    '''
    # put coordinates in correct order
    coordinates = np.array(coords).T

    # create figure and display
    fig=plt.figure()
    ax=fig.add_subplot(projection='3d')
    ax.plot3D(coordinates[0],coordinates[1],coordinates[2],ltype)

    return
        


###############################################################################
def distanceMatrix( xyz ):
    '''
    use broadcasting to calculate map efficiently. in following, indices
    are numbered from (1) fastest to slowest.
    the shape of the N residue xyz array is A=(N(2),3(1))
    the shape of the higher rank created array is B=(N(3),:(2),3(1))
    then DMV = B-A. explanation: B's index(2) is stretched to N, 
    to B=(N(3),N(2),3(1)). and then A is expanded to rank 3 and stretched to
    A=(N(3), N(2), 3(1)). Note that if we consider this a matrix of 3-vectors
    that is NxN, we have A has one column per coordinate (all same along
    the new stretched slow 3-axis), and B has one row per coordinate (all same
    along the new stretched slower 2-axis).

    Args:
        xyz (array): coordinates, shape = (N,3)

    Returns:
        (array, float) distance matrix, shape = (N,N)

    '''                                         
    DMV = xyz[:,np.newaxis,:] - xyz  
    
    # now sum squares along fast axis and take square root
    return np.sqrt( (DMV*DMV).sum(axis=2) )  


###############################################################################
def contactMatrix( xyz, cutoff=5.0 ):
    '''
    see distanceMatrix() for description of distance calculation.
    contacts are determined using cut-off.

    Args:
        xyz (array, float): coordinates, shape = (N,3)
        cutoff (float, optional): cut off distance. Defaults to 5.0.
        
    Returns:
        (array, bool) contact matrix, shape = (N,N)

    '''                                         
    DMV = xyz[:,np.newaxis,:] - xyz  
    DM = np.sqrt( (DMV*DMV).sum(axis=2) ) 
    
    # now create the contact ap
    return np.less(DM,cutoff)

###############################################################################
def heatMap( squareMatrix, title='Heat Map', label='position',cmap='afmhot_r',\
            log=False, tickLabels=None):

    heatMapDict = {'variable': squareMatrix,
                   'name': label,
                   'title': title,
                   'colormap': cmap
                   }
    if log: heatMapDict['variable'] = np.log(heatMapDict['variable'])

    fig,ax=plt.subplots()
    ax.imshow(heatMapDict['variable'])
    ax.set_xlabel( heatMapDict['name'] +' 1' )
    ax.set_ylabel( heatMapDict['name'] +' 2' )
    ax.xaxis.set_label_position("top")
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    
    if tickLabels:
        ax.set_xticks(range(20),tickLabels, rotation='vertical')
        ax.set_yticks(range(20),tickLabels, rotation='horizontal')
       
    ax.imshow( heatMapDict['variable'], cmap=heatMapDict['colormap'])

    return





###############################################################################
def aaContactMatrix( xyz, seq, cutoff=5.0 ):
    '''
    uses a singular similarity transform to transform the sequence based 
    contact matrix (see contactMatrix for details) to an amino-acid type
    based matrix. The 'sequence matrix' ( shape=(n,20) ) has a single
    non-zero element (=1) in each row (column k for amino-acid type k).
    note 1: the uncorrected transformed contact matrix will count self-contacts
    once, and all other contacts twice. For type-k to type-k contacts, both counts
    add to diagonal, where as for k to l!=k contacts, they are split on the 
    symmetric positions (+1 at k,l and +1 at l,k).
    the algorithm corrects the calculated sequence based contact matrix 
    (subtracting the identity) to remove self contacts before transformation.
    note 2: double counting still results in double countacts along the
    diagonal----should be corrected ...?
    
    Args:
        xyz (array, float): coordinates, shape = (N,3)
        seq (array, float or int): the encoded sequence, shape = (N)
        cutoff (float, optional): cut off distance. Defaults to 5.0.

    Returns:
        TYPE: DESCRIPTION.

    '''
    # calculate distance matrix
    DMV = xyz[:,np.newaxis,:] - xyz  
    DM = np.sqrt( (DMV*DMV).sum(axis=2) )  

    # calculate contact matrix & remove self-contacts
    CM = np.less(DM,cutoff)
    np.fill_diagonal(CM, False)

    # calculate sequence matrix
    seqMatrix = np.zeros( ( len(seq), 20 ) )
    seqMatrix[ range(len(seq)), seq] = 1

    # perform similarity transform on contact map to xform to aa based map
    return np.linalg.multi_dot([seqMatrix.T,CM,seqMatrix])




