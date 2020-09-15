#!/usr/bin/env python
# coding: utf-8

prot_id = "1eej.pdb"
prot_file = "1eej.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]

residuescationA = []
residuesPIA = []
listeA=[]

parser=argparse.ArgumentParser(
    description='''Ce programme identifie les liaisons cations-pi à partir d'un fichier Protein Data Bank (PDB). Les critères pris en compte proviennent du Protein Interaction Calculator que l'on peut retrouver en suivant le lien : http://pic.mbu.iisc.ernet.in/PIC_Criteria.pdf. Le parser de Biopython est strucutré de la manière suivante : Structure/model/chain/residu/atome Fonctions utilisées :''',
    epilog="""Je vous souhaite une bonne utilisation du programme.""")
parser.add_argument('parser.get_structure', type=int, default=42, help='Creation of a structure object from a PDB file')
parser.add_argument('objet.get_name', type=int, default=43, help='Renvoie le nom correspondant à l objet : Structure/model/chain/residu/atome')
parser.add_argument('atom.get_id', type=int, default=44, help='Renvoie le numéro rattaché au résidue dans le fichier PDB')
args=parser.parse_args()

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
                    v=(r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    listeA.append(v)
                    print("possible cation-pi interactions in chain A:",r1.get_parent(),r2.get_parent()),
                    print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    print (round(r1[r1name] - r2[r2name],2))
                
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
                    v=(r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    listeA.append(v)
                    print("possible cation-pi interactions in chain :", r1.get_parent(),r2.get_parent()),
                    print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    print (round(r1[r1name] - r2[r2name],2)) 

                

h=0                   
while h< len(listeA):
    a=listeA[h]
    x=listeA.count(a)
    if x>=3:
        print(a)
    h=h+1

