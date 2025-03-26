#!\usr\bin\env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 07:59:14 2024

@author: allen

Preprocessing of mmCIF file to extract DNA chains and store coordinate
data in array format designed for easy use in structure analysis.

"""
import os
import numpy as  np
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

structureDirectory = '/home/allen/projects/DATA/db/assemblies'
structureFile = '1qrv-assembly1.cif'
outputDirectory = 'data'
logFile = 'mSA.log'


'''
###############################################################################
########################  FUNCTIONS  ######################################
###############################################################################
'''
phosphateAtoms = { 'P', 'OP1', 'OP2', 'OP3' }
riboseAtoms = { "C1\'", "C2\'", "C3\'", "C4\'", "C5\'", "O3\'", "O4\'", "O5\'"}
baseAtoms = { 'N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'N2', \
             'O6', 'N6' }
dnaResidues = {'DA', 'DC','DG', 'DT'}
backboneAtoms={'N','CA','O','OXT'} # for now, include CA (for glycine)
nucleotides = "ACGT"
aminoAcids = "ARNDCEQGHILKMFPSTWYV"

def extractDnaStructure( chain ):
    
    xyzPhosphate=[]; xyzRibose=[]; xyzBase=[]; seq=[]
    for residue in chain:
        if residue.get_resname() not in dnaResidues: continue # skip if not DNA   
        
        #######################################################################
        phArray=np.ones((4,3))*np.nan # arrays to receive substructure coords
        rbArray=np.ones((8,3))*np.nan
        baArray=np.ones((20,3))*np.nan
        ph=[]; rb=[]; ba=[]  # lists for substructure coords
        for atom in residue:
            if atom.get_name() in phosphateAtoms: ph.append(atom.get_coord())
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
        
    return xyzPhosphate, xyzRibose, xyzBase, seq

def extractProteinStructure( chain ):
    
    xyzBackbone=[]; xyzSidechain=[]; seq=[]
    for residue in chain:
        if not is_aa(residue): continue   # skip if not aa residue
             
        #######################################################################
        # arrays to receive substructure coords
        bbArray=np.ones((4,3))*np.nan; scArray=np.ones((11,3))*np.nan
        bb=[]; sc=[] # lists for substructures
        for atom in residue:
            if atom.get_name() in backboneAtoms: bb.append(atom.get_coord())
            else: sc.append(atom.get_coord())
                
        # fill arrays with coordinate lists
        bbArray[:len(bb)] = bb; scArray[:len(sc)] = sc
        #######################################################################

        #######################################################################
        # create the oneHot array for the sequence
        seqArray = np.zeros(20)
        index=aminoacids.find(seq1(residue.get_resname())) # use 1-letter code
        seqArray[index] = 1
        #######################################################################
        
        # add arrays to chain lists
        xyzBackbone.append(bbArray)
        xyzSidechain.append(scArray)
        seq.append(seqArray)
        
    return xyzBackbone, xyzSidechain, seq 
        
 
    
    
    
    


'''
###############################################################################
########################  MAIN  ######################################
###############################################################################
'''
# N.B. MMCIFParser assigned author chainids
# the sequence determined is from the residues present in the structure
# the seq of crystallized protein can be gotten from:
#_entity_poly.pdbx_seq_one_letter_code_can 

logOut = open(logFile,'w')
'''
1. load and parse the mmCIF file into a biopythn structure object
'''
# assume pdb id code given by first 4 characters of file name
code = structureFile[:4]
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure(code,os.path.join(structureDirectory,structureFile))

# for log file
# list number of models, and number of chains in model 0
logOut.write(f'{code}: {len(structure)} model(s). ')
chains = []
model = structure[0]
logOut.write(\
   f'{len(model)} chain(s) in model 0: {tuple(model.child_dict.keys())} \n')

'''
2. setup dictionaries for protein and dna data
'''
# create the structure and sequence arrays
# each entry in the dictionary will be a structure/sequence for a chain
phosphateDict = {}; riboseDict={}; baseDict = {}; dnaSequenceDict = {}
backboneDict = {}; sidechainDict={}; proteinSequenceDict = {}

for chain in model:   # loop through all chains in model 0
    logOut.write(f'processing structure chain {chain.id} ')
    
    if chain is dna:   
        
        phosphate, ribose, base, seq = extractDnaStructure( chain )
        if phosphate and ribose and base:   # only add if all non-empty
            phosphateDict[chain.id] = np.array(phosphate)
            riboseDict[chain.id] = np.array(ribose)
            baseDict[chain.id] = np.array(base)
            logOut.write('structure added. ')
        
        if seq:
            dnaSequenceDict[chain.id] = np.array(seq)
            logOut.write('sequence added.\n')
                           
    elif chain is protein:
        logOut.write('as protein. ')
        xyzBackbone, xyzSidechain, seq = extractProteinStructure( chain )   
        if xyzBackbone and xyzSidechain:   # only add if all non-empty
            backboneDict[chain.id] = np.array(xyzBackbone)
            sidechainDict[chain.id] = np.array(xyzSidechain)
            logOut.write('structure added. ')
        if seq:
            proteinSequenceDict[chain.id] = np.array(seq)
            logOut.write('sequence added.\n')
        

    '''
os.makedirs(outputDirectory,exist_ok=True)
for k in baseDict.keys():
    filePath = os.path.join(outputDirectory,code+'_'+k)
    print('saving structure chain', k,':',filePath)
    np.savez(filePath,seq=oneHotDict[k], ph=phosphateDict[k],\
             rb=riboseDict[k], ba=baseDict[k] )
             
    '''         
logOut.close()
