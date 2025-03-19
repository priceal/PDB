#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:50:36 2024

@author: allen
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import MMCIFParser, is_aa
from Bio.SeqUtils import seq1
import os
import numpy as  np

structureDirectory = '../DATA/db/assemblies'
structureFile = '1a66-assembly1.cif'

###############################################################################


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

ca=[]
# Iterate of all chains in the model in order to find all residues
for chain in model:
  # Iterate of all residues in each model in order to find proper atoms
  for res in chain:
    # Check if residue number ( .get_id() ) is in the list
    if is_aa(res):
      # Append CA atom to list
      ca.append(res['CA'].get_coord())
      
      



ca=np.array(ca)

cas = np.array(ca).T

fig=plt.figure()
ax=fig.add_subplot(projection='3d')
ax.plot3D(cas[0],cas[1],cas[2],'-')
