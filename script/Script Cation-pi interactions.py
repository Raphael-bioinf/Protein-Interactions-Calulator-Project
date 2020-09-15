#!/usr/bin/env python
# coding: utf-8

# In[28]:


#PHYDROPHOBIC INTERACTIONS:The following residues are considered to participate in interactions if they fall within 5Å range.
#ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, TYR.

prot_id = "1eej.pdb"
prot_file = "1eej.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]

residuescationA = []
residuesPIA = []
residuescationB = []
residuesPIB = []
listeA=[]
listeB=[]

x=0
z=0
for chain in model:
    a=""
    b=""
    c="chain"
    cationAA=["LYS","ARG"]
    PIAA=["PHE","TRP","TYR"]
    if x==0 :
        lettres = chr(ord('a')+x)
        a=model[lettres.capitalize()]
        for residue in a:
            if residue.get_resname() in cationAA:
                residuescationA.append(residue)
            if residue.get_resname() in PIAA:
                residuesPIA.append(residue)
    
    if x==1:
        lettres = chr(ord('a')+x)
        b=model[lettres.capitalize()]
        
        for residue in b:
            if residue.get_resname() in cationAA:
                residuescationB.append(residue)
            if residue.get_resname() in PIAA:
                residuesPIB.append(residue)
    x=x+1

    
for r1 in residuescationA:
    for atom in r1:
        r1name=atom.get_name()
        r1index = r1.get_id()[1]
    for r2 in residuesPIA:
        for atom in r2:
            r2name=atom.get_name()
            r2index = r2.get_id()[1]
            if(r1index != r2index and len(r1name)>1 and r1name != "CA" and len(r2name)>1 and r2name != "CA" and "H" not in r1name and "H" not in r2name):
                try:
                    distance = r1[r1name] - r2[r2name]
                except KeyError:
                    continue
                if (r1[r1name] - r2[r2name])< 6:
                    v=(r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    listeA.append(v)
                    #print("possible cation-pi interactions in chain A:"),
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
                try:
                    distance = r1[r1name] - r2[r2name]
                except KeyError:
                    continue
                if (r1[r1name] - r2[r2name])< 6:
                    v=(r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    listeA.append(v)
                    #print("possible cation-pi interactions in chain A:"),
                    #print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                    #print (round(r1[r1name] - r2[r2name],2))            
                
for r3 in residuescationB:
    for atom in r3:
        r3name=atom.get_name()
        r3index = r3.get_id()[1]
    for r4 in residuesPIB:
        for atom in r4:
            r4name=atom.get_name()
            r4index = r4.get_id()[1]
            if(r3index != r4index and len(r3name)>1 and r3name != "CA" and len(r4name)>1 and r4name != "CA" and "H" not in r3name and "H" not in r4name):
                try:
                    distance = r3[r3name] - r4[r4name]
                except KeyError:
                    continue
                if (r3[r3name] - r4[r4name])< 6:
                    u=(r3.get_resname(),r3index,"-",r4.get_resname(),r4index,)
                    listeB.append(u)
                    #print("possible cation-pi interactions in chain B:"),
                    #rint (r3.get_resname(),r3index,"-",r4.get_resname(),r4index,)
                    #print (round(r3[r3name] - r4[r4name],2))
for r3 in residuesPIB:
    for atom in r3:
        r3name=atom.get_name()
        r3index = r3.get_id()[1]
    for r4 in residuescationB:
        for atom in r4:
            r4name=atom.get_name()
            r4index = r4.get_id()[1]
            if(r3index != r4index and len(r3name)>1 and r3name != "CA" and len(r4name)>1 and r4name != "CA" and "H" not in r3name and "H" not in r4name):
                try:
                    distance = r3[r3name] - r4[r4name]
                except KeyError:
                    continue
                if (r3[r3name] - r4[r4name])< 6:
                    u=(r3.get_resname(),r3index,"-",r4.get_resname(),r4index,)
                    listeB.append(u)
                   # print("possible cation-pi interactions in chain B:"),
                    #print (r3.get_resname(),r3index,"-",r4.get_resname(),r4index,)
                    #print (round(r3[r3name] - r4[r4name],2))
h=0                   
while h< len(listeA):
    a=listeA[h]
    x=listeA.count(a)
    if x>=4:
        print(a)
    h=h+1
while h< len(listeB):
    a=listeB[h]
    x=listeB.count(a)
    if x>=4:
        print(a)
    h=h+1


# In[ ]:




