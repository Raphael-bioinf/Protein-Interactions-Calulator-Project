#!/usr/bin/env python
# coding: utf-8

# In[1]:


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]


hydrogen = 
reidues = []
for chain in model:
    for res in chain:
        if res.get_resname() in hydrogen:
            residues.append(res)

    
def hbondfun(acceptor, donnor, distON = 3.5, distS = 4):
    d = acceptor - donnor
    if "O" in acceptor.get_name() or "N" in acceptor.get_name():
        if d < distON:
            return(d)
    elif "S" in acceptor.get_name():
        if d < distS:
            return(d)


def hbond_main_mainfun(atom1, atom2, distON = 3.5, distS = 4):
    name1 = atom1.get_name()
    name2 = atom2.get_name()
    if len(name1) == 1 and len(name2) == 1:
        if name1 == "O" and name2 == "H":
            d = hbondfun(
                acceptor=atom1,
                donnor=atom2,
                distON=distON,
                distS=distS)
            if d:
                return(d)
        elif name2 == "O" and name1 == "H":
            d = hbondfun(
                acceptor = atom2,
                donnor = atom1,
                distON = distON,
                distS = distS)
            if d:
                return(d)


def hbond_main_sidefun(atom1, atom2, distON = 3.5, distS = 4):
    name1 = atom1.get_name()
    name2 = atom2.get_name()
    if len(name1) == 1 and len(name2) == 2:
        if (name1 in ["O", "N", "S"]) and (
                "O" in name2 or "N" in name2 or "S" in name2):
            d = hbondfun(
                acceptor = atom1,
                donnor = atom2,
                distON = distON,
                distS = distS)
            if d:
                return(d)
    if len(name1) == 2 and len(name2) == 1:
        if (name2 in ["O", "N", "S"]) and (
                "O" in name1 or "N" in name1 or "S" in name1):
            d = hbondfun(
                acceptor = atom1,
                donnor = atom2,
                distON = distON,
                distS = distS)
            if d:
                return(d)


def hbond_side_sidefun(atom1, atom2, distON = 3.5, distS = 4):
    name1 = atom1.get_name()
    name2 = atom2.get_name()
    if len(name1) == 2 and len(name2) == 2:
        if ("O" in name2 or "N" in name2 or "S" in name2) and (
                "O" in name1 or "N" in name1):
            d = hbondfun(
                acceptor = atom2,
                donnor = atom1,
                distON = distON,
                distS = distS)
            if d:
                return(d)
        elif ("O" in name1 or "N" in name1 or "S" in name1) and ("O" in name2 or "N" in name2):
            d = hbondfun(
                acceptor = atom1,
                donnor = atom2,
                distON = distON,
                distS = distS)
            if d:
                return(d)







