#!/usr/bin/env python
# coding: utf-8

# In[1]:


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

from Bio.PDB import Atom
from math import sqrt
from Bio.PDB import PDBParser

parser = PDBParser(PERMISSIVE = 1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]


arosul = ["PHE", "TRP",  "TYR", "CYS", "MET"] # AA impliqués dans des interactions cycle aromatique - soufre
residues = []


for chain in model:
    for res in chain:
        if res.get_resname() in arosul:
            residues.append(res)


def dist_cal(x1, y1, z1, x2, y2, z2): # fonction calculant la distance euclidienne entre 2 points de coordonnées x, y & z
    dist = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    return(dist)


def center_mass(resid): # calcul des coordonnées des centres des cycles
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
        resid.center_mass = (xmean, ymean, zmean)
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


dist = 5.3
for res1 in residues:
    for res2 in residues:
        sulres = ["CYS", "MET"]
        
        if (res1.get_resname() in sulres and res2.get_resname() not in sulres):
            for atom in res1:
                if "S" in atom.get_name():
                    s_coord = atom.get_coord()
                    d = dist_cal(
                        s_coord[0],
                        s_coord[1],
                        s_coord[2],
                        center_mass(res2).center_mass[0], 
                        center_mass(res2).center_mass[1], 
                        center_mass(res2).center_mass[2]) 
                    if d < dist:
                        print(res1.get_resname(), res1.get_id()[1], res2.get_resname(), res2.get_id()[1], "dist =", d)
        elif (res2.get_resname() in sulres and res1.get_resname() not in sulres):
            for atom in res2:
                if "S" in atom.get_name():
                    s_coord = atom.get_coord()
                    d = dist_cal(
                        s_coord[0],
                        s_coord[1],
                        s_coord[2],
                        center_mass(res1).center_mass[0], 
                        center_mass(res1).center_mass[1], 
                        center_mass(res1).center_mass[2])
                    if d < dist:
                        print(res1.get_resname(), res1.get_id()[1], res2.get_resname(), res2.get_id()[1], "dist =", d)