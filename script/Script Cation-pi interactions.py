#!/usr/bin/env python
# coding: utf-8

# In[21]:

import sys

prot_id = "4rcn.pdb"
prot_file = sys.argv[1] #"4rcn.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]

if "-h" in sys.argv or "--help" in sys.argv:
    print ("Ce programme identifie les liaisons cations-pi à partir d'un fichier Protein Data Bank (PDB). Les critères pris en compte proviennent du Protein Interaction Calculator que l'on peut retrouver en suivant le lien : http://pic.mbu.iisc.ernet.in/PIC_Criteria.pdf. Le parser de Biopython est strucutré de la manière suivante : Structure/model/chain/residu/atome.")
    print ("Fonctions utilisées :")
    print ('parser.get_structure -->  ', 'Creation of a structure object from a PDB file')
    print ('objet.get_name -->  ', 'Renvoie le nom correspondant à l objet : Structure/model/chain/residu/atome')
    print ('parser.get_structur -->e  ', 'Renvoie le numéro rattaché au résidue dans le fichier PDB')
    print ('')
    
residuescationA = []
residuesPIA = []
listeA=[]


for chain in model:
    cationAA=["LYS","ARG"]
    PIAA=["PHE","TRP","TYR"]
    for residue in chain:
        if residue.get_resname() in cationAA:
                residuescationA.append(residue)
        if residue.get_resname() in PIAA:
                residuesPIA.append(residue)
   
    
for r1 in residuescationA:
    for atom in r1:
        r1name=atom.get_name()
        r1index = r1.get_id()[1]
    for r2 in residuesPIA:
        for atom in r2:
            r2name=atom.get_name()
            r2index = r2.get_id()[1]
            if(r1index != r2index and len(r1name)>1 and r1name != "CA" and len(r2name)>1 and r2name != "CA" and "H" not in r1name and "H" not in r2name):
                if (r1[r1name] - r2[r2name])< 6:
                    v=(r1.get_resname(),r1index,r1.get_parent(),"-",r2.get_resname(),r2index,r2.get_parent())
                    listeA.append(v)
                    #print("possible cation-pi interactions in chain A:",r1.get_parent(),r2.get_parent()),
                    #print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    #print (round(r1[r1name] - r2[r2name],2))
                
for r1 in residuesPIA:
    for atom in r1:
        r1name=atom.get_name()
        r1index = r1.get_id()[1]
    for r2 in residuescationA:
        for atom in r2:
            r2name=atom.get_name()
            r2index = r2.get_id()[1]
            if(r1index != r2index and len(r1name)>1 and r1name != "CA" and len(r2name)>1 and r2name != "CA" and "H" not in r1name and "H" not in r2name):
                if (r1[r1name] - r2[r2name])< 6:
                    v=(r1.get_resname(),r1index,r1.get_parent(),"-",r2.get_resname(),r2index,r2.get_parent())
                    listeA.append(v)
                    #print("possible cation-pi interactions in chain :", r1.get_parent(),r2.get_parent()),
                    #print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    #print (round(r1[r1name] - r2[r2name],2)) 

                

h=0                   
while h< len(listeA):
    a=listeA[h]
    x=listeA.count(a)
    if x>=5:
        print(a)
    h=h+1







