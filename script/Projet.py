#!/usr/bin/env python
# coding: utf-8
# In[]:

# Programme permettant d'accéder aux données et aux coordonnées du fichier PDB


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]
for chain in model:
    for residue in chain:
        for atom in residue:
            print (chain, residue, atom, atom.get_coord())







