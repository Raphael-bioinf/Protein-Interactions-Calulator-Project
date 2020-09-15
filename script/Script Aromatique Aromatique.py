#!/usr/bin/env python
# coding: utf-8

# In[1]:


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

from math import sqrt
from Bio.PDB import Atom
from Bio.PDB import PDBParser

parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]

achain = model['A']
bchain = model['B']
    

aroaro = ["PHE",  "TRP", "TYR"]
residues = []
for res in achain:
    if res.get_resname() in aroaro:
        residues.append(res)
#print(residues)


def dist_cal(x1, y1, z1, x2, y2, z2): # fonction calculant la distance euclidienne entre 2 points de coordonnées x, y & z
    dist = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    return(dist)


def dist_center_mass_calc(resid1, resid2):
    xmean = 0
    ymean = 0
    zmean = 0
    
    if resid1.get_resname() in ["PHE", "TYR"]:
        for atom in resid1:
            if atom.get_name() in ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]:
                xmean = xmean + atom.get_coord()[0]
                ymean = ymean + atom.get_coord()[1]
                zmean = zmean + atom.get_coord()[2]
        xmean = xmean / 6
        ymean = ymean / 6
        zmean = zmean / 6
        resid1_center_mass = (xmean, ymean, zmean)
    elif resid1.get_resname() == "TRP":
        for atom in resid1:
            if atom.get_name() in ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]:
                xmean = xmean + atom.get_coord()[0]
                ymean = ymean + atom.get_coord()[1]
                zmean = zmean + atom.get_coord()[2]
        xmean = xmean / 6
        ymean = ymean / 6
        zmean = zmean / 6
    resid1_center_mass = (xmean, ymean, zmean)
    
    if resid2.get_resname() in ["PHE", "TYR"]:
        for atom in resid2:
            if atom.get_name() in ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]:
                xmean = xmean + atom.get_coord()[0]
                ymean = ymean + atom.get_coord()[1]
                zmean = zmean + atom.get_coord()[2]
        xmean = xmean / 6
        ymean = ymean / 6
        zmean = zmean / 6
        resid2_center_mass = (xmean, ymean, zmean)
    elif resid2.get_resname() == "TRP":
        for atom in resid2:
            if atom.get_name() in ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]:
                xmean = xmean + atom.get_coord()[0]
                ymean = ymean + atom.get_coord()[1]
                zmean = zmean + atom.get_coord()[2]
        xmean = xmean / 6
        ymean = ymean / 6
        zmean = zmean / 6
    resid2_center_mass = (xmean, ymean, zmean)
    
    dist = dist_cal(
        x1 = resid1_center_mass[0],
        y1 = resid1_center_mass[1],
        z1 = resid1_center_mass[2],
        x2 = resid2_center_mass[0],
        y2 = resid2_center_mass[1],
        z2 = resid2_center_mass[2])
    return(dist)


def aroarofun(resid1, resid2, dmin = 4.5, dmax = 7): # détermine si une distance correspond à une interaction entre AA aromatique

    d = dist_center_mass_calc(resid1, resid2)

    if (d > dmin) and (d < dmax):
        print(resid1.get_resname(), resid1.get_id()[1], resid2.get_resname(), resid2.get_id()[2], "dist =", d)
        return(d)


for res1 in residues:
    for res2 in residues:
        aroarofun(res1, res2)