#!/usr/bin/env python
# coding: utf-8

import sys

prot_id = "4rcn.pdb"
prot_file = sys.argv[1] #"4rcn.pdb"

#On va utiliser le parser de Biopython qui nous permet d'accéder aux éléments d'un fichier PDB.
from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]

#On définit la fonction help qui peut être appelé pour une présentation du programme et des fonctions utilisées
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

#Cette première boucle définit et filtre les acides aminés d'intérêts pour les intéractions cation-pi
for chain in model:
    cationAA=["LYS","ARG"]
    PIAA=["PHE","TRP","TYR"]
    for residue in chain:
        if residue.get_resname() in cationAA:
                residuescationA.append(residue)
        if residue.get_resname() in PIAA:
                residuesPIA.append(residue)
   
#Cette boucle nous permet de sélectionner les atomes responsables des interactions cations-pi puis de filtrer ceux respectant notre critère de distance.
#Ils sont ensuite stocké dans la une liste
for r1 in residuescationA:
    cation=['NH1','NH2','NZ']
    for atom in r1:
        r1name=atom.get_name()
        if r1name in cation:
            r1index = r1.get_id()[1]
    for r2 in residuesPIA:
        for atom in r2:
            r2name=atom.get_name()
            r2index = r2.get_id()[1]
            if(r1index != r2index and len(r1name)>1 and r1name != "CA" and r1name != "CB"):
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
        cation=['NH1','NH2','NZ']
        for atom in r2:
            r2name=atom.get_name()
            if r2name in cation:
                r2index = r2.get_id()[1]
            if(r1index != r2index and len(r1name)>1 and r1name != "CA" and r1name != "CB"):
                if (r1[r1name] - r2[r2name])< 6:
                    v=(r1.get_resname(),r1index,r1.get_parent(),"-",r2.get_resname(),r2index,r2.get_parent())
                    listeA.append(v)
                    #print("possible cation-pi interactions in chain :", r1.get_parent(),r2.get_parent()),
                    #print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    #print (round(r1[r1name] - r2[r2name],2)) 

                
#cette boucle while nous permet de ne ressortir que les résidus présentant 6 interactions différentes ( pour les 6 carbones du cycle aromatique)
h=0                   
while h< len(listeA):
    a=listeA[h]
    x=listeA.count(a)
    if x>=6:
        print(a)
    h=h+1







