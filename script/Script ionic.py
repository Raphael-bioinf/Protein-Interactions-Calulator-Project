#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio.PDB import PDBParser
import sys


prot_id = "5AGY.pdb"
prot_file = sys.argv[1]


parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]


if "-h" in sys.argv or "--help" in sys.argv:
    print("Ce programme identifie les interactions ioniques à partir d'un fichier Protein Data Bank (PDB). Les critères pris en compte proviennent du Protein Interaction Calculator que l'on peut retrouver en suivant le lien : http://pic.mbu.iisc.ernet.in/PIC_Criteria.pdf. Le parser de Biopython est strucutré de la manière suivante : Structure/model/chain/residu/atome.")
    print("Fonctions utilisées :")
    print('parser.get_structure -->  ', 'Creation of a structure object from a PDB file')
    print('objet.get_name -->  ', 'Renvoie le nom correspondant à l objet : Structure/model/chain/residu/atome')
    print('parser.get_structure -->  ', 'Renvoie le numéro rattaché au résidue dans le fichier PDB')
    print('')

    
ionic = ["ARG",  "LYS", "HIS", "ASP", "GLU"]
residues = []


for chain in model: # protéine -> chaîne -> résidues impliqués dans interactions ioniques
    for res in chain:
        if res.get_resname() in ionic:
            residues.append(res)


def ionic_fun(atom1, atom2, resid1, resid2, dist = 6):
    if (
        (("N" in atom1.get_name() and len(atom1.get_name()) > 1
          and resid1.get_resname() in ["LYS", "ARG", "HIS"])
         and ("O" in atom2.get_name() and len(atom2.get_name()) > 1
              and resid2.get_resname() in ["ASP", "GLU"])) or
        (("O" in atom1.get_name() and len(atom1.get_name()) > 1
          and resid1.get_resname() in ["ASP", "GLU"])
         and ("N" in atom2.get_name() and len(atom2.get_name()) > 1
              and resid2.get_resname() in ["LYS", "ARG", "HIS"]))
    ):
        d = atom1 - atom2
        if d < dist:
            print(resid1.get_resname(), resid1.get_id()[1],resid2.get_resname(), resid2.get_id()[1], "dist =", d)
            return(d)


for res1 in residues:
    for res2 in residues:
        for atom1 in res1:
            for atom2 in res2:
                ionic_fun(atom1, atom2, res1, res2, 6)