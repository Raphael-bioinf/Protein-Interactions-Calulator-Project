#!/usr/bin/env python
# coding: utf-8

# In[1]:

import sys

prot_id = "5AGY.pdb"
prot_file = sys.argv[1]

from Bio.PDB import PDBParser
parser = PDBParser(PERMISSIVE = 1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]


if "-h" in sys.argv or "--help" in sys.argv:
    print("Ce programme identifie les liaisons hydrogènes à partir d'un fichier Protein Data Bank (PDB). Les critères pris en compte proviennent du Protein Interaction Calculator que l'on peut retrouver en suivant le lien : http://pic.mbu.iisc.ernet.in/PIC_Criteria.pdf. Le parser de Biopython est strucutré de la manière suivante : Structure/model/chain/residu/atome.")
    print("Fonctions utilisées :")
    print('parser.get_structure -->  ', 'Creation of a structure object from a PDB file')
    print('objet.get_name -->  ', 'Renvoie le nom correspondant à l objet : Structure/model/chain/residu/atome')
    print('parser.get_structure -->  ', 'Renvoie le numéro rattaché au résidue dans le fichier PDB')
    print('')


residues = []
for chain in model:
    for res in chain:
        residues.append(res)

    
def hbond_fun(acceptor, donnor, resid1, resid2, distON = 3.5, distS = 4):
    d = acceptor - donnor
    if resid1.get_id()[1] != resid2.get_id()[1] and resid1.get_resname() != resid2.get_resname():
        if "O" in donnor.get_name() or "N" in donnor.get_name() or "S" in donnor.get_name():
            if "O" in acceptor.get_name() or "N" in acceptor.get_name():
                if d < distON:
                    print(resid1.get_id()[1], resid1.get_resname(), acceptor.get_name(), resid2.get_id()[1], resid2.get_resname(), donnor.get_name(), "dist =", d)
                    return(d)
    elif resid1.get_id()[1] == resid2.get_id()[1]:
        pass


for res1 in residues:
    for res2 in residues:
        for acceptor in res1:
            for donnor in res2:
                hbond_fun(acceptor, donnor, res1, res2)