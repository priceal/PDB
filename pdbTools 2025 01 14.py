#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:46:56 2024

@author: allen
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import requests


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
def fetch( code ):
    '''
    downloads PDB file and returns a list of strings, one entry
    per line

    Args:
        code (str): the 4 character PDB code

    Returns:
        lines (list): the text of the file.

    '''
    url = 'https://files.rcsb.org/download/'+code+'.pdb'
    download = requests.get(url)
    lines = download.text.split('\n')
    
    return lines

###############################################################################
def load( filename ):
    '''
    reads in a PDB file, returns list same as fetch

    Args:
        filename (str): path to file

    Returns:
        lines (list): the text of the file.

    '''
    with open(filename,'r') as file:
        lines = file.readlines()
        
    return lines
        
###############################################################################
def getCoordinates( source, chain='A', select='all'):
    '''
    extracts coordinates from source
    
    Args:
        source (str or list): path to the PDB file or a list containing
            the lines of the PDB file
        chain (str, optional): which chain to read. Defaults to 'A'.
        select (string, optional): choices are 'all' or 'ca'. 
            Defaults to 'all'.

    Returns:
        numpy array of coordinates, shape = (N,3)

    '''    
    # if filename given, read in file and store each line in pdbData
    if type(source) is str:
        with open(source,'r') as file:
            lines = file.readlines()
    else: lines=source   # assume list of pdb file given
    
   # extract atom coordinates
    coords = []
    if select=='ca':
        for line in lines:
            lineSplit = line.split()
            if len(lineSplit)>0 and lineSplit[0] == 'ATOM' and \
                lineSplit[4]==chain and lineSplit[2]=='CA':
                coords.append([float(x) for x in lineSplit[6:9]])   
'''
    elif select=='sc':
        residue=0; com=np.zeros(3); atomCount=0
        for line in lines:
            lineSplit = line.split()
            if len(lineSplit)>0 and lineSplit[0] == 'ATOM' and \
                lineSplit[4]==chain:
                newResidue = int(lineSplit[5])
                if newResidue==residue or residue==0:
                    com += np.array(([float(x) for x in lineSplit[6:9]])) 
                    atomCount += 1
                    residue=newResidue
                else:
                    coords.append(com/atomCount)
                    com=np.zeros(3)
                    atomCount=0
 '''                   
   
    elif select=='all':
        for line in lines:
            lineSplit = line.split()
            if len(lineSplit)>0 and lineSplit[0] == 'ATOM' and \
                lineSplit[4]==chain:
                coords.append([float(x) for x in lineSplit[6:9]])
   
    return np.array(coords)

##############################################################################
def getSequence( source, resCode = defaultCode, chain='A' ):
    '''
    extracts sequence and returns coded version

    Args:
        source (str or list): path to the PDB file or a list containing
            the lines of the PDB file
        resCode (dict, optional): which endoding to use.
            Defaults to defaultCode.
        chain (str, optional): which chain to read. Defaults to 'A'.

    Returns:
        (array, int) encoded sequence, shape = (N)

    '''
    # if filename given, read in file and store each line in pdbData
    if type(source) is str:
        with open(source,'r') as file:
            pdbData = file.readlines()
    else: pdbData=source   # assume list of pdb file given

    # extract sequence code
    sequence = []
    for line in pdbData:
        lineSplit = line.split()
        if len(lineSplit)>0 and lineSplit[0] == 'ATOM' and \
            lineSplit[4]==chain and lineSplit[2]=='CA':
            sequence.append( resCode[ lineSplit[3] ] )  
     
    return np.array(sequence)

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
