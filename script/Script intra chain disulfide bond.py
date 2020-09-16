# Problème à régler : exprime en duplicata les interactions (A avec B et B avec A)

import sys

prot_id = "1eej.pdb"
prot_file = sys.argv[1]

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]

cysresiduesA = []
cysresiduesB = []
x=0
for chain in model:
    a=""
    b=""
    c="chain"
    if x==0 :
        lettres = chr(ord('a')+x)
        a=model[lettres.capitalize()]
        for residue in a:
            if residue.get_resname() == "CYS":
                cysresiduesA.append(residue)
    
    if x==1:
        lettres = chr(ord('a')+x)
        b=model[lettres.capitalize()]
        
        for residue in b:
            if residue.get_resname() == "CYS":
                cysresiduesB.append(residue)
    x=x+1
                
    
for c1 in cysresiduesA:
    c1index = c1.get_id()[1]
    for c2 in cysresiduesA:
        c2index = c2.get_id()[1]
        if(c1index != c2index):
            if(c1['SG'] - c2['SG'])< 2.2:
                print("possible di-sulfide bond in chain A:"),
                print ("Cys",c1index,"-")
                print ("Cys",c2index,)
                print (round(c1['SG'] - c2['SG'],2))
for c3 in cysresiduesB:
    c3index = c3.get_id()[1]
    for c4 in cysresiduesB:
        c4index = c4.get_id()[1]
        if(c3index != c4index):
            if (c3['SG'] - c4['SG'])< 2.2:
                print("possible di-sulfide bond in chain B:"),
                print ("Cys",c3index,"-")
                print ("Cys",c4index,)
                print (round(c3['SG'] - c4['SG'],2))







