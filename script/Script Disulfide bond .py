#!/usr/bin/env python
# coding: utf-8

# In[1]:

import sys

prot_id = "5AGY.pdb"
prot_file = sys.argv[1]

#On va utiliser le parser de Biopython qui nous permet d'accéder aux éléments d'un fichier PDB.
from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]

#Le parser de Biopython est strucutré de manière Structure/model/chain/residu/atome
cysresidues = []
for chain in model:
    for residue in chain:
        # On sélectionne uniquement les résidue de Cystéine responsable des liaison disulfide
        if residue.get_resname() == "CYS":
            cysresidues.append(residue)
    for c1 in cysresidues:
        #La fonction get_id retourne le numéro affecté àl'acide aminé dans le fichier PDB
        c1index = c1.get_id()[1] 
        for c2 in cysresidues:
            c2index = c2.get_id()[1]
            if(c1index!=c2index):
                 # On conditionne la présence d'une liaison disulfide à une distance de 2.2 Angström
                if(c1['SG'] - c2['SG'])< 2.2:
                    print("possible di-sulfide bond:",r1.get_parent(),r2.get_parent())
                    print ("Cys",c1index,"-")
                    print ("Cys",c2index,)
                    print (round(c1['SG'] - c2['SG'],2))







