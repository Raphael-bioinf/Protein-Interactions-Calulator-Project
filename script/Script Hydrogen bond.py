#!/usr/bin/env python
# coding: utf-8

# In[1]:


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

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


