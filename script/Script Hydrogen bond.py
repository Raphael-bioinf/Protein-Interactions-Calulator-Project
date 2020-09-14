#!/usr/bin/env python
# coding: utf-8

# In[1]:


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]
achain = model['A']
bchain = model['B']
    
elecgiver = []
for residue in achain:
    if residue.get_resname() == "CYS":
        elecgiver.append(residue)
for c1 in elecgiver:
    c1index = c1.get_id()[1]
    for c2 in elecgiver:
        c2index = c2.get_id()[1]
        if(c1['SG'] - c2['SG'])< 2.2:
            print("possible di-sulfide bond:"),
            print ("Cys",c1index,"-")
            print ("Cys",c2index,)
            print (round(c1['SG'] - c2['SG'],2))







