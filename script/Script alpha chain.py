#!/usr/bin/env python
# coding: utf-8




prot_id = "5AGY.pdb"
prot_file = "5AGY.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]
achain = model['A']
bchain = model['B']

for residue in achain:
    index = residue.get_id()[1]
    calpha = residue['CA']
    carbon = residue['C']
    nitrogen = residue['N']
    oxygen = residue['O']
    print ("residue:", residue.get_resname(),index)
    print ("N - CA", (nitrogen - calpha))
    print ("CA - C", (calpha - carbon))
    print ("C - O", (carbon - oxygen))
    
    if achain.has_id(index+1):
        nextresidue = achain[index+1]
        nextnitrogen = nextresidue['N']
        print ("C - N", (carbon - nextnitrogen))
    print()







