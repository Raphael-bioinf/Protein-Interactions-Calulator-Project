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
#ALA(),    
#ILE(,,
#LEU(,, ),
#MET(, CE, CG,SD),  
#PHE(,  , 
#PRO(,    
#TRP(,     
#TYR(,  
#VAL(,

                                         
for r1 in residuesA:
    cresidues =['CA','CB','CD1','CD2','CG1','CG2','CE1','CE2','CG','CD2','CG','CE3','CH2','CZ2','CZ3','CZ','CE']
    r1index = r1.get_id()[1]
    for r2 in residuesA:
        r2index = r2.get_id()[1]
        if(r1index != r2index):
            if(r1['CB'] - r2['CB'])< 5:
                print("possible hydrophobic interactions in chain A:"),
                print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                print (round(r1['CB'] - r2['CB'],2))
            try:
                distance = r1['CD1'] - r2['CD1']
            except KeyError:
                    continue
            if(r1['CD1'] - r2['CD1'])< 5:
                print("possible hydrophobic interactions in chain A:"),
                print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                print (round(r1['CD1'] - r2['CD1'],2))
            try:
                distance = r1['CD2'] - r2['CD2']
            except KeyError:
                    continue
            if(r1['CD2'] - r2['CD2'])< 5:
                print("possible hydrophobic interactions in chain A:"),
                print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                print (round(r1['CD2'] - r2['CD2'],2))
                
for r3 in residuesB:
    r3index = r3.get_id()[1]
    for r4 in residuesB:
        r4index = r4.get_id()[1]
        if(r3index != r4index):
            if (r3['CB'] - r4['CB'])< 5:
                print("possible hydrophobic interactionsin chain B:"),
                print (r3.get_resname(),r3index,"-",r4.get_resname(),r4index,)
                print (round(r3['CB'] - r4['CB'],2))
            try:
                distance = r3['CD1'] - r4['CD1']
            except KeyError:
                    continue
            if (r3['CD1'] - r4['CD1'])< 5:
                print("possible hydrophobic interactionsin chain B:"),
                print (r3.get_resname(),r3index,"-",r4.get_resname(),r4index,)
                print (round(r3['CD1'] - r4['CD1'],2))
            try:
                distance = r3['CD2'] - r4['CD2']
            except KeyError:
                    continue
            if (r3['CD2'] - r4['CD2'])< 5:
                print("possible hydrophobic interactionsin chain B:"),
                print (r3.get_resname(),r3index,"-",r4.get_resname(),r4index,)
                print (round(r3['CD2'] - r4['CD2'],2))

