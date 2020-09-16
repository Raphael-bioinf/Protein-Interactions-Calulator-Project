#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio.PDB import PDBParser
import sys


prot_id = "5AGY.pdb"
prot_file = sys.argv[1]


parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]

    
ionic = ["ARG",  "LYS", "HIS", "ASP", "GLU"]
residues = []


for chain in model: # protéine -> chaîne -> résidues impliqués dans interactions ioniques
    for res in chain:
        if res.get_resname() in ionic:
            residues.append(res)


def center_mass(resid):
    xmean = 0
    ymean = 0
    zmean = 0
    if resid.get_resname() in ["PHE", "TYR"]:
        for atom in resid:
            if atom.get_name() in ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]:
                xmean = xmean + atom.get_coord()[0]
                ymean = ymean + atom.get_coord()[1]
                zmean = zmean + atom.get_coord()[2]
        xmean = xmean / 6
        ymean = ymean / 6
        zmean = zmean / 6
    elif resid.get_resname() == "TRP":
        for atom in resid:
            if atom.get_name() in ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]:
                xmean = xmean + atom.get_coord()[0]
                ymean = ymean + atom.get_coord()[1]
                zmean = zmean + atom.get_coord()[2]
        xmean = xmean / 6
        ymean = ymean / 6
        zmean = zmean / 6
    resid.center_mass = (xmean, ymean, zmean)
    return(resid)


def ionicfun(atom1, atom2, resid1, resid2, dist = 6):
    if (
        (("N" in atom1.get_name() and len(atom1.get_name()) > 1
          and resid1.get_resname() in ["LYS", "ARG", "HIS"])
         and ("O" in atom2.get_name() and len(atom2.get_name()) > 1
              and resid2.get_resname() in ["ASP", "GLU"])) or
        (("O" in atom1.get_name() and len(atom1.get_name()) > 1
          and resid1.get_resname() in ["ASP", "GLU"])
         and ("N" in atom2.get_name() and len(atom2.get_name()) > 1
              and resid2.get_resname() in ["LYS", "ARG", "HIS"]))
    ):
        d = atom1 - atom2
        if d < dist+0.3 and d > dist-0.3:
            print(resid1.get_resname(), resid1.get_id()[1],resid2.get_resname(), resid2.get_id()[1], "dist =", d)
            return(d)


for res1 in residues:
    for res2 in residues:
        for atom1 in res1:
            for atom2 in res2:
                ionicfun(atom1, atom2, res1, res2, 6)