# -*- coding: utf-8 -*-
"""
Calculates and plots as a heat map the distance matrix for two sets of
structures. input files must be .npz files, with keys that match those in 
inputs below.

coordinate arrays must have shape (N,maxResSize,3) where N is number of 
residues, maxResSize is number of atoms/res. Therefore, indices indicate 
(residue #, atom #, coordinate).
"""
import os
import numpy as np
import matplotlib.pyplot as plt

'''
###############################################################################
######################### functions ###########################################
###############################################################################
'''

def loadStructureAndSequence( files, group, Directory='' ):
    
    listCoords = []  # to receive all coords from list of files
    listSeq = []  # to receive sequence from files
    for f in files:   
        if Directory:                   # use directory if given
            f=os.path.join(Directory,f) 
            
        tempStructure=np.load(f)    
        coords = []  # to receive coords from groups from one file
        for g in group:
            coords.append(tempStructure[g])
            
        # concatenate along residue atom number axis (1) and append
        listCoords.append(np.concatenate(coords,axis=1))
        listSeq.append(tempStructure['seq']) 
        
    return np.concatenate(listCoords,axis=0), \
           np.concatenate(listSeq,axis=0)
           
###############################################################################

def distanceMatrix( structure1, structure2 ):
    
    '''
    first calculate the distance map. use broadcasting to calculate map efficiently. 
    the shape of structure1 array is S1=(N1,R1,3) and that of structure2 is
    S2=(N2,R2,3). New axes are created to increase the ranks to that of
    S1*=(N1,1,R1,1,3) and S2*=(1,N2,1,R2,3). Broadcasting will stretch the axes 
    of dim=1 thus: S1* -> (N1,(N2),R1,(R2),3) and S2* -> ((N1),N2,(R1),R2,3).
    The shape of the resulting difference will be 
    S1*-S2* = (N1(0),N2(1),R1(2),R2(3),3(4)). where the indices meanings are:
        
        (0) residue from structure1
        (1) residue from structure2
        (2) atom from structure1 residue
        (3) atom from structure2 residue
        (4) coordinate x, y, or z
    
    To get distance between atoms, take the sqrt of sums along axis 4.
    To get minimum distance between residues,take minimum along axes 2 and 3.
    '''

    # the vector difference between the two structures
    vDM=structure1[          :, np.newaxis,           :, np.newaxis, : ] - \
        structure2[ np.newaxis,           :, np.newaxis,          :, : ]
      
    # the atom-atom distance matrix
    aDM = np.sqrt( (vDM*vDM).sum( axis = 4 ) )  # 4 is axis of coordinates (xyz)
    
    # the min-distance res-res matrix
    return np.nanmin( aDM, axis=(2,3))    # (2,3) are axes of atom # w/in residues
   
###############################################################################

def displayHeatMap( heatMapDict, xticks='', yticks='' ):
    
    if heatMapDict[ 'logarithmic' ]==True:
        heatMapDict[ 'variable' ] = np.log(heatMapDict[ 'variable' ])
            
    fig,ax=plt.subplots()
    ax.imshow(heatMapDict['variable'])
    ax.set_xlabel( heatMapDict['xlabel'] )
    ax.set_ylabel( heatMapDict['ylabel'] )
    ax.xaxis.set_label_position("top")
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    if heatMapDict['title']:
        ax.set_title( heatMapDict['title'] )
    else:
        ax.set_title('heat map')
        
    if xticks:
        ax.set_xticks(range(len(xticks)))
        ax.set_xticklabels(xticks)
    if yticks:
        ax.set_yticks(range(len(yticks)))
        ax.set_yticklabels(yticks)
        
    ax.imshow( heatMapDict['variable'], cmap=heatMapDict['colormap'])
    return

'''
###############################################################################
############################# main ############################################
###############################################################################
'''

if __name__ == "__main__":
    '''
    structure file information --- format of dictionary entries:
        label : [ [files], {groups} ]
    
    label = used to label axes
    [files] = list of .npz files for this set of structure
    {groups} = set of subgroups of atoms to consider
    allowed groups ...
               proteins:   bb = backbone, sc = sidechain
               dna:        ph = phosphate, rb = ribose, ba = base
    
    N.B. all .npz files under same label MUST have structure arrays of same shape!           
    '''

    inputStructures = { 'protein':  [ ['1le5_A.npz','1le5_B.npz'],    \
                                     {'sc'} ],
                        'dna':      [ ['1le5_C.npz','1le5_D.npz'],    \
                                     {'ba'} ]
                      }
 
    # optional structure file directory, can leave undefined '' or '.'
    fileDirectory = '../DATA/PDB/npz'
    
    cutoff = 0  # non-zero for a contact map with cutoff value
    mapTitle = 'sidechain-base contacts'
    colorMap = 'OrRd'
    logorithmic = True    
    
    ###############################################################################
    
    structure = {}; sequence = {}; polymerType = {}    
    for k,(files,group) in inputStructures.items():
        if group&{'sc','bb'}:
            polymerType[k] = 'protein'
        elif group&{'ph','rb','ba'}:
            polymerType[k] = 'dna'    
        structure[k], sequence[k] = loadStructureAndSequence(
                                    files, group, 
                                    Directory=fileDirectory
                                    )
    # the vector displacement matrix: shape = (N1,N2,R1,R2,3)
    label2, structure2 = structure.popitem() # pull out coordinates LIFO
    label1, structure1 = structure.popitem()
    DM = distanceMatrix(structure1, structure2)
    
    # now create the contact map if asked
    if cutoff: DM = np.less(DM,cutoff)
    
    displayHeatMap( {'variable': DM,
                     'title': mapTitle,
                     'xlabel': label2,
                     'ylabel': label1,
                     'colormap': colorMap,
                     'logarithmic': logorithmic
                     } 
                  )
    
