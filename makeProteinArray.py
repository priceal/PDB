#!\usr\bin\env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 07:59:14 2024

@author: allen
"""
from Bio.PDB import MMCIFParser, is_aa
from Bio.SeqUtils import seq1
import os
import numpy as  np

structureDirectory = '/home/allen/projects/DATA/db/assemblies'
structureFile = '1b72-assembly1.cif'
code = '1b72'

outputDirectory = 'data'

###############################################################################
# N.B. MMCIFParser assigned author chainids
# the sequence determined is from the residues present in the structure
# the seq of crystallized protein can be gotten from:
#_entity_poly.pdbx_seq_one_letter_code_can 

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
backboneDict = {}; sidechainDict={}; sequenceDict = {}
backboneAtoms={'N','CA','O','OXT'} # for now, do not include CA (for glycine)
for chainid in chains:   # loop through all chains in model 0
    print('processing structure chain', chainid)
    xyzBackbone=[]; xyzSidechain=[]; seq=[]
    for residue in model[chainid]:
        if not is_aa(residue):   # skip if not aa residue
            continue
        
        bbArray=np.ones((4,3))*np.nan; scArray=np.ones((11,3))*np.nan
        bb=[]; sc=[]
        for atom in residue:
            if atom.get_name() in backboneAtoms: bb.append(atom.get_coord())
            else: sc.append(atom.get_coord())
                
        bbArray[:len(bb)] = bb; scArray[:len(sc)] = sc
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

# create the oneHot array for the sequence
aminoAcids = "ARNDCEQGHILKMFPSTWYV"
oneHotDict = {}
for k, v in sequenceDict.items():
    print('chain',k)
    pssm=np.zeros((len(v),20))
    pssm[ range(len(v)), [aminoAcids.find(c) for c in v]  ] =1
    oneHotDict[k] = pssm
print('completed processing all chains')

os.makedirs(outputDirectory,exist_ok=True)
for k in sidechainDict.keys():
    filePath = os.path.join(outputDirectory,code+'_'+k)
    print('saving structure chain', k,':',filePath)
    np.savez(filePath,seq=oneHotDict[k],bb=backboneDict[k],sc=sidechainDict[k])

