#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#PHYDROPHOBIC INTERACTIONS:The following residues are considered to participate in interactions if they fall within 5Å range.
#ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, TYR.

prot_id = "1eej.pdb"
prot_file = "1eej.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]

residuesA = []
residuesB = []
x=0
for chain in model:
    a=""
    b=""
    c="chain"
    hydrophobicAA=["ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","TYR"]
    if x==0 :
        lettres = chr(ord('a')+x)
        a=model[lettres.capitalize()]
        for residue in a:
            if residue.get_resname() in hydrophobicAA:
                residuesA.append(residue)
    
    if x==1:
        lettres = chr(ord('a')+x)
        b=model[lettres.capitalize()]
        
        for residue in b:
            if residue.get_resname() in hydrophobicAA:
                residuesB.append(residue)
    x=x+1

#For protein backbones, this representation preserves five backbone atoms for each amino acid: nitrogen (N), the alpha carbon (CA), the carbonyl carbon (C), the carbonyl oxygen (O), and the polar hydrogen on nitrogen. The sidechain is replaced by the beta carbon (CB) 
#and a CEN atom whose radius and properties (polarity, charge, etc.) are determined by the residue's identity
cresidues =['CA','CB','O','N','C'] 
for r1 in residuesA:
    for atom in r1:
        r1name=atom.get_name()
        if r1name in cresidues:
            r1index = r1.get_id()[1]
    for r2 in residuesA:
        for atom in r2:
            r2name=atom.get_name()
            if r2name in cresidues:
                r2index = r2.get_id()[1]
            if(r1index != r2index):
                try:
                    distance = r1[r1name] - r2[r2name]
                except KeyError:
                    continue
                if 3.4<(r1[r1name] - r2[r2name])< 5:
                    print("possible hydrophobic interactions in chain A:"),
                    print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    print (round(r1[r1name] - r2[r2name],2))
            
                
for r3 in residuesB:
    for atom in r3:
        r3name=atom.get_name()
        if r3name in cresidues:
            r3index = r3.get_id()[1]
    for r4 in residuesB:
        for atom in r4:
            r4name=atom.get_name()
            if r4name in cresidues:
                r4index = r4.get_id()[1]
            if(r3index != r4index):
                try:
                    distance = r3[r3name] - r4[r4name]
                except KeyError:
                    continue
                if 3.4<(r3[r3name] - r4[r4name])< 5:
                    print("possible hydrophobic interactions in chain B:"),
                    print (r3.get_resname(),r3index,"-",r4.get_resname(),r4index,)
                    print (round(r3[r3name] - r4[r4name],2))
           

