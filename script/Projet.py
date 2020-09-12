#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import scipy as sp

prot_id = "5AGY.pdb"
prot_file = "5AGY.pdb"
arr_x = [];
arr_y = [];

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
for model in structure.get_list():
    for chain in model.get_list():
        for residue in chain.get_list():
            for atom in residue:
                y = residue.get_atoms()
                y = list(y)
                arr_y.append({'Atoms': [y[0]]})
                
for model in structure.get_list():
    for chain in model.get_list():
        for residue in chain.get_list():
            for atom in residue:
                x = atom.get_coord()
                arr_x.append({'X': [x[0]],'Y':[x[1]],'Z':[x[2]]})
  
        
sample = pd.DataFrame(arr_x)
sample2 = pd.DataFrame(arr_y)
sample2[''] = arr_x

print(sample2)  



  


# In[4]:


import pandas as pd
import numpy as np
import scipy as sp

prot_id = "5AGY.pdb"
prot_file = "5AGY.pdb"
arr_x = [];
arr_y = [];
arr_z = [];
from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
for model in structure.get_list():
    for chain in model.get_list():
        for residue in chain.get_list():
            for atom in model:
                y = model.get_atoms()
                y = list(y)
                arr_y.append({'Atoms': [y[0]]})
for model in structure.get_list():
    for chain in model.get_list():
        for residue in chain.get_list():
            for atom in residue:
                x = atom.get_coord()
                arr_x.append({'X': [x[0]],'Y':[x[1]],'Z':[x[2]]})
  
        
sample = pd.DataFrame(arr_x)
sample2 = pd.DataFrame(arr_y)


print(sample2)  


# In[8]:


prot_id = "5AGY.pdb"
prot_file = "5AGY.pdb"
arr_x = [];
arr_y = [];

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]
for chain in model:
    for residue in chain:
        for atom in residue:
            print (chain, residue, atom, atom.get_coord())


# In[ ]:




