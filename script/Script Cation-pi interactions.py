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

