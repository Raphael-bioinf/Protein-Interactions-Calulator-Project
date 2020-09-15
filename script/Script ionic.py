#!/usr/bin/env python
# coding: utf-8

# In[1]:


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

from Bio.PDB import PDBParser
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]

achain = model['A']
#bchain = model['B']
    
ionic = ["ARG",  "LYS", "HIS", "ASP", "GLU"]
residues = []
for res in achain:
    if res.get_resname() in ionic:
        residues.append(res)


def dist_cal(x1, y1, z1, x2, y2, z2):
    dist = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    return(dist)


def dist_center_mass_calc(resid1, resid2):
    dist = dist_cal(
        x1=resid1.center_mass[0],
        y1=resid1.center_mass[1],
        z1=resid1.center_mass[2],
        x2=resid2.center_mass[0],
        y2=resid2.center_mass[1],
        z2=resid2.center_mass[2])
    return(dist)


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
        if d < dist:
            print(resid1.get_resname(), resid1.get_id()[1],resid2.get_resname(), resid2.get_id()[1], "dist =", d)
            return(d)


reverse_residues = residues[::-1]
for res1 in residues:
    for res2 in reverse_residues:
        for atom1 in res1:
            for atom2 in res2:
                ionicfun(atom1, atom2, res1, res2)