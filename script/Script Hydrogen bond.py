#!/usr/bin/env python
# coding: utf-8

# In[1]:


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

from Bio.PDB import PDBParser
parser = PDBParser(PERMISSIVE = 1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]


hydrogen = []
residues = []
for chain in model:
    for res in chain:
        residues.append(res)

    
def hbondfun(acceptor, donnor, distON = 3.5, distS = 4):
    d = acceptor - donnor
    if "O" in acceptor.get_name() or "N" in acceptor.get_name():
        if d < distON:
            print(acceptor.get_name(), donnor.get_name(), "dist =", d)
            return(d)
    elif "S" in acceptor.get_name():
        if d < distS:
            print(acceptor.get_name(), donnor.get_name(), "dist =", d)
            return(d)


for res in residues:
    for atom1 in res:
        for atom2 in res:
            hbondfun(atom1, atom2)


