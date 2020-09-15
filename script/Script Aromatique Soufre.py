#!/usr/bin/env python
# coding: utf-8

# In[1]:


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"

from Bio.PDB import Atom
from Bio.PDB import PDBParser

parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]

achain = model['A']
bchain = model['B']


arosul = ["PHE", "TRP",  "TYR", "CYS", "MET"]
residues = []
for res in achain:
    if res.get_resname() in arosul:
        residues.append(res)
#print(residues)


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
                        res2.center_mass[0], #need to correct "AttributeError: 'Residue' object has no attribute 'center_mass'" error
                        res2.center_mass[1], #need to correct "AttributeError: 'Residue' object has no attribute 'center_mass'" error
                        res2.center_mass[2]) #need to correct "AttributeError: 'Residue' object has no attribute 'center_mass'" error
                    if d < dist:
                        print(d)
        elif (res2.get_resname() in sulres and res1.get_resname() not in sulres):
            for atom in res2:
                if "S" in atom.get_name():
                    s_coord = atom.get_coord()
                    d = dist_cal(
                        s_coord[0],
                        s_coord[1],
                        s_coord[2],
                        res1.center_mass[0], #need to correct "AttributeError: 'Residue' object has no attribute 'center_mass'" error
                        res1.center_mass[1], #need to correct "AttributeError: 'Residue' object has no attribute 'center_mass'" error
                        res1.center_mass[2]) #need to correct "AttributeError: 'Residue' object has no attribute 'center_mass'" error
                    if d < dist:
                        print(d)