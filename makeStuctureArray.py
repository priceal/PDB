#!\usr\bin\env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 07:59:14 2024

@author: allen

Preprocessing of mmCIF file to extract DNA chains and store coordinate
data in array format designed for easy use in structure analysis.

"""
from Bio.PDB import MMCIFParser
import os
import numpy as  np

structureDirectory = '/home/allen/projects/DATA/db/assemblies'
structureFile = '1qrv-assembly1.cif'
outputDirectory = 'data'



phosphateAtoms = { 'P', 'OP1', 'OP2' }
riboseAtoms = { "C1\'", "C2\'", "C3\'", "C4\'", "C5\'", "O3\'", "O4\'", "O5\'"}
baseAtoms = { 'N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'N2', \
             'O6', 'N6' }
dnaResidues = {'DA', 'DC','DG', 'DT'}
backboneAtoms={'N','CA','O','OXT'} # for now, do not include CA (for glycine)


def extractDnaStructure( chain ):
    
    xyzPhosphate=[]; xyzRibose=[]; xyzBase=[]; seq=[]
    for residue in chain:
        if residue.get_resname() not in dnaResidues: continue # skip if not DNA       
        phArray=np.ones((3,3))*np.nan # arrays to receive substructure coords
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
        
        # now add substructure arrays to chain lists
        xyzPhosphate.append(phArray) 
        xyzRibose.append(rbArray)
        xyzBase.append(baArray)
        seq.append(residue.get_resname()[1:])
        
        
    
def extractProteinStructure( chain ):
    
    xyzBackbone=[]; xyzSidechain=[]; seq=[]
    for residue in chain:
        if not is_aa(residue):   # skip if not aa residue
            continue
        
        # arrays to receive substructure coords
        bbArray=np.ones((4,3))*np.nan; scArray=np.ones((11,3))*np.nan
        bb=[]; sc=[] # lists for substructures
        for atom in residue:
            if atom.get_name() in backboneAtoms: bb.append(atom.get_coord())
            else: sc.append(atom.get_coord())
                
        # fill arrays with coordinate lists
        bbArray[:len(bb)] = bb; scArray[:len(sc)] = sc
        
        # add arrays to chain lists
        xyzBackbone.append(bbArray); xyzSidechain.append(scArray)
        seq.append(residue.get_resname())
        
    if xyzBackbone and xyzSidechain:   # only add to dict if an list is non-empty
        backboneDict[chainid] = np.array(xyzBackbone)
        sidechainDict[chainid] = np.array(xyzSidechain)
        print('    added protein structure', end=', ')
        
    if seq:
        sequenceDict[chainid] = seq1(''.join(seq))
        print('    added protein sequence')
        
print('completed processing all chains')









###############################################################################
# N.B. MMCIFParser assigned author chainids
# the sequence determined is from the residues present in the structure
# the seq of crystallized protein can be gotten from:
#_entity_poly.pdbx_seq_one_letter_code_can 

# assume pdb id code given by first 4 characters of file name
code = structureFile[:4]

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure(code,os.path.join(structureDirectory,structureFile))

# list number of models, and number of chains in model 0
print(len(structure),'model(s)')
chains = []
model = structure[0]
for chain in model:
    chains.append( chain.get_id() )
print(len(model),'chain(s) in model 0:',chains)

# create the protein structure and sequence arrays
# each entry in the dictionary will be a structure/sequence for a chain
phosphateDict = {}; riboseDict={}; baseDict = {}
sequenceDict = {}


# create the protein structure and sequence arrays
# each entry in the dictionary will be a structure/sequence for a chain
backboneDict = {}; sidechainDict={}; sequenceDict = {}

for chainid in chains:   # loop through all chains in model 0
    print('processing structure chain', chainid)
    
    
    
    
    if chain is dna:
        
        xyzPhosphate, xyzRibose, xyzBase, seq = \
            extractDnaStructure( chain )
            
            
    if xyzPhosphate and xyzRibose and xyzBase:   # only add if all non-empty
        phosphateDict[chainid] = np.array(xyzPhosphate)
        riboseDict[chainid] = np.array(xyzRibose)
        baseDict[chainid] = np.array(xyzBase)

        print('    added DNA structure', end=', ')
        
    if seq:
        sequenceDict[chainid] = ''.join(seq)
        print('    added DNA sequence')
                
        print('completed processing all chains')

        # create the oneHot array for the sequence
        nucleotides = "ACGT"
        oneHotDict = {}
        for k, v in sequenceDict.items():
            print('chain',k)
            pssm=np.zeros((len(v),4))
            pssm[ range(len(v)), [nucleotides.find(c) for c in v]  ] =1
            oneHotDict[k] = pssm
        print('completed processing all chains')

            
            
    elif chain is protein:
        
        xyzBackbone, xyzSidechain, seq=[] = \
            extractProteinStructure( chain )   
            
       
    
    
os.makedirs(outputDirectory,exist_ok=True)
for k in baseDict.keys():
    filePath = os.path.join(outputDirectory,code+'_'+k)
    print('saving structure chain', k,':',filePath)
    np.savez(filePath,seq=oneHotDict[k], ph=phosphateDict[k],\
             rb=riboseDict[k], ba=baseDict[k] )
             
             
