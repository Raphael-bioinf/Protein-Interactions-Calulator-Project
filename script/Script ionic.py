#!/usr/bin/env python
# coding: utf-8

# In[1]:


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

from Bio.PDB import PDBParser
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]

achain = model['A']
bchain = model['B']
    
ionicres = []
for residue in achain:
    if (residue.get_resname() == "ARG" or residue.get_resname() == "LYS" or residue.get_resname() == "HIS" or residue.get_resname() == "ASP" or residue.get_resname() == "GLU"):
        ionicres.append(residue)
        #print(ionicres)
        
        for c1 in ionicres:
            c1index = c1.get_id()[1]
            
            for c2 in ionicres:
                c2index = c2.get_id()[1]
                
                if(c1['OE2'] - c2['N'] < 8 and c1['OE2'] - c2['N'] > 4):
                    print("possible ionic bond:"),
                    print (residue,c1index,"-")
                    print (residue,c2index,)
                    print (round(c1['OE2'] - c2['N'],2))