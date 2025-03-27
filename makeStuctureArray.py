#!\usr\bin\env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 07:59:14 2024

@author: allen

Preprocessing of mmCIF files to extract protein/DNA chains and store 
coordinate data in array format designed for easy use in structure analysis.

"""
import os
import numpy as  np
import pandas as pd
from Bio.PDB import MMCIFParser, is_aa
from Bio.SeqUtils import seq1

'''
###############################################################################
########################   FUNCTIONS   ########################################
###############################################################################
'''

# variable definitions
unexceptedAtoms = { 'H', 'ZN' }
phosphateAtoms = { 'P', 'OP1', 'OP2', 'OP3' }
riboseAtoms = { "C1\'", "C2\'", "C3\'", "C4\'", "C5\'", "O3\'", "O4\'", "O5\'"}
baseAtoms = { 'N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'N2', \
             'O6', 'N6' }
dnaResidues = {'DA', 'DC','DG', 'DT'}
backboneAtoms={'N','CA','O','OXT'} # for now, include CA (for glycine)
nucleotides = "ACGT"
aminoAcids = "ARNDCEQGHILKMFPSTWYV"

###############################################################################
def is_dna( chain ):
    '''
    determines if at least 2 residues in chain are canonical DNA residue types

    Args:
        chain (TYPE): DESCRIPTION.

    Returns:
        bool: DESCRIPTION.

    '''
    count = 0
    for r in chain:
        if r.resname in dnaResidues:
            count += 1
            
    return count > 1

###############################################################################
def is_protein( chain ):
    '''
    determines if at least 5 residues in chain are canonical AA residue types 

    Args:
        chain (TYPE): DESCRIPTION.

    Returns:
        bool: DESCRIPTION.

    '''
    count = 0
    for r in chain:
        count += is_aa(r)
            
    return count > 4

###############################################################################
def extractDnaStructure( chain ):
    '''
    create structure and sequence arrays from a protein chain

    Parameters
    ----------
    chain : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    xyzPhosphate=[]; xyzRibose=[]; xyzBase=[]; seq=[]
    for residue in chain:
        if residue.get_resname() not in dnaResidues: continue # skip if not DNA   
        
        #######################################################################
        phArray=np.ones((4,3))*np.nan # arrays to receive substructure coords
        rbArray=np.ones((8,3))*np.nan
        baArray=np.ones((20,3))*np.nan
        ph=[]; rb=[]; ba=[]  # lists for substructure coords
        for atom in residue:
            if atom.element in unexceptedAtoms: pass
            elif atom.get_name() in phosphateAtoms: ph.append(atom.get_coord())
            elif atom.get_name() in riboseAtoms: rb.append(atom.get_coord())
            elif atom.get_name() in baseAtoms: ba.append(atom.get_coord())
            else: pass  # usually this is for H in NMR stuctures
            # issue: non-canonical bases with specil atom types
        
        # fill arrays with coordinate lists up to length of list, force
        # correct shape of array conversion of lists
        phArray[:len(ph)] = np.array(ph).reshape(-1,3)
        rbArray[:len(rb)] = np.array(rb).reshape(-1,3)
        baArray[:len(ba)] = np.array(ba).reshape(-1,3)    
        #######################################################################
        
        #######################################################################
        # create the oneHot array for the sequence
        seqArray = np.zeros(4)
        index=nucleotides.find(residue.get_resname()[1:]) # drop the 'D'
        seqArray[index] = 1
        #######################################################################
       
        # now append substructure arrays to chain lists of residues
        xyzPhosphate.append(phArray) 
        xyzRibose.append(rbArray)
        xyzBase.append(baArray)
        seq.append(seqArray)
        
    return np.array(xyzPhosphate), np.array(xyzRibose), \
            np.array(xyzBase), np.array(seq)

###############################################################################
def extractProteinStructure( chain ):
    '''
    create structure and sequence arrays from a DNA chain   

    Parameters
    ----------
    chain : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    xyzBackbone=[]; xyzSidechain=[]; seq=[]
    for residue in chain:
        if not seq1(residue.get_resname()): continue # skip if empty
        if seq1(residue.get_resname()) not in aminoAcids: continue # skip non-AA 
             
        #######################################################################
        # arrays to receive substructure coords
        bbArray=np.ones((4,3))*np.nan; scArray=np.ones((11,3))*np.nan
        bb=[]; sc=[] # lists for substructures
        for atom in residue:
            if atom.element in unexceptedAtoms: pass
            elif atom.get_name() in backboneAtoms: bb.append(atom.get_coord())
            else: sc.append(atom.get_coord())
                
        # fill arrays with coordinate lists
        try:
            bbArray[:len(bb)] = bb
            scArray[:len(sc)] = sc
        except:
            print(chain,residue,residue.get_resname(),end=' ')
            for a in residue:
                print(a.element,end=' ')
            print('\n')
        #######################################################################

        #######################################################################
        # create the oneHot array for the sequence
        seqArray = np.zeros(20)
        index=aminoAcids.find(seq1(residue.get_resname())) # use 1-letter code
        seqArray[index] = 1
        #######################################################################
        
        # add arrays to chain lists
        xyzBackbone.append(bbArray)
        xyzSidechain.append(scArray)
        seq.append(seqArray)
        
    return np.array(xyzBackbone), np.array(xyzSidechain), np.array(seq) 

'''
###############################################################################
########################  MAIN  ######################################
###############################################################################
'''

# variables for main
# inputs
#structureDirectory = '/home/allen/projects/DATA/db/assemblies'
structureListFile = './csv/sort715.csv'

#outputs
outputDirectory = '/home/allen/projects/DATA/PDB/npz'
logFile = 'mSA_715.log'

# N.B. MMCIFParser assigns author chainids
# the sequence determined is from the residues present in the structure
# the seq of crystallized protein can be gotten from:
#_entity_poly.pdbx_seq_one_letter_code_can 

# initialize I/O
logOut = open(logFile,'w')
structureDf = pd.read_csv( structureListFile )
logOut.write(f'structure list loaded: {structureListFile}\n')
parser = MMCIFParser(QUIET=True)

# main loop
for entry in structureDf.itertuples():
    print(entry.pdbid,end=' ')

    '''
    1. load and parse the mmCIF file into a biopython structure object
    '''
    
    structure = parser.get_structure(entry.pdbid,entry.path)
    
    # for log file
    # list number of models, and number of chains in model 0
    logOut.write(f'\nfile: {entry.path} \n')
    logOut.write(f'{entry.pdbid}: {len(structure)} model(s). ')
    chains = []
    model = structure[0]
    logOut.write(\
       f'{len(model)} chain(s) in model 0: {tuple(model.child_dict.keys())} \n')
    
    '''
    2. make data directory if needed
    '''
    # create the structure and sequence arrays
    # each entry in the dictionary will be a structure/sequence for a chain
    os.makedirs(outputDirectory,exist_ok=True)
    
    '''
    3. loop through all chains and send to correct extractor
    '''
    dnaChains = []
    proteinChains = []
    for chain in model:   # loop through all chains in model 0
        logOut.write(f'processing chain {chain.id} ')
        filePath = os.path.join(outputDirectory,entry.pdbid+'_'+chain.id)
        
        if is_dna( chain ):   
            phosphate, ribose, base, seq = extractDnaStructure( chain )
            logOut.write('dna     ')
            np.savez(filePath,seq=seq, ph=phosphate, rb=ribose, ba=base )
            logOut.write(f': {filePath}.npz\n')           
            dnaChains.append(chain.id)
                               
        elif is_protein( chain ):
            backbone, sidechain, seq = extractProteinStructure( chain )  
            logOut.write('protein ')
            np.savez(filePath,seq=seq,bb=backbone,sc=sidechain)
            logOut.write(f': {filePath}.npz\n')
            proteinChains.append(chain.id)
            
        else:
            logOut.write('***TYPE NOT IDENTIFIED***\n')
    
    # add classification of chains to structureDf
    structureDf.at[entry.Index,'dna chains'] = ','.join(dnaChains)
    structureDf.at[entry.Index,'protein chains'] = ','.join(proteinChains)

# resave new dataframe w/ updated columns listing chain classifications
structureDf.to_csv(structureListFile,index=False)                   
logOut.close()
