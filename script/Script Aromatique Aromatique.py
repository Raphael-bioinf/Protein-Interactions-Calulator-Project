#!/usr/bin/env python
# coding: utf-8

from Bio.PDB import Atom
from math import sqrt
from Bio.PDB import PDBParser
import argparse
import sys


prot_id = "5AGY.pdb"
prot_file = sys.argv[1]


parser = PDBParser(PERMISSIVE = 1)
structure = parser.get_structure(prot_id, prot_file)
model = structure[0]


if "-h" in sys.argv or "--help" in sys.argv:
    print("Ce programme identifie les interactions entre cycles aromatiques à partir d'un fichier Protein Data Bank (PDB). Les critères pris en compte proviennent du Protein Interaction Calculator que l'on peut retrouver en suivant le lien : http://pic.mbu.iisc.ernet.in/PIC_Criteria.pdf. Le parser de Biopython est strucutré de la manière suivante : Structure/model/chain/residu/atome.")
    print("Fonctions utilisées :")
    print('parser.get_structure -->  ', 'Creation of a structure object from a PDB file')
    print('objet.get_name -->  ', 'Renvoie le nom correspondant à l objet : Structure/model/chain/residu/atome')
    print('parser.get_structure -->  ', 'Renvoie le numéro rattaché au résidue dans le fichier PDB')
    print('')


residues = []
aroaro = ["PHE", "TRP", "TYR"]
for chain in model: # protéine -> chaîne -> résidues impliqués dans interactions aromatiques / aromatiques
    for res in chain:
        if res.get_resname() in aroaro:
            residues.append(res)


def dist_cal(x1, y1, z1, x2, y2, z2): # fonction calculant la distance euclidienne entre 2 points de coordonnées x, y & z
    dist = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    return(dist)


def dist_center_coord(resid1, resid2): # calcul coordonnées centre cycles aromatiques
    xmean = 0
    ymean = 0
    zmean = 0
    
    if resid1.get_resname() in ["PHE", "TYR"]: # cycle 1
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
    
    if resid2.get_resname() in ["PHE", "TYR"]: # cycle 2
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


def aro_aro(resid1, resid2, dmin = 4.5, dmax = 7): # détermine si une distance correspond à une interaction entre AA aromatique

    d = dist_center_coord(resid1, resid2)

    if (d > dmin) and (d < dmax):
        print( resid1.get_id()[1],resid1.get_resname(), resid1.get_parent(), resid2.get_id()[1], resid2.get_resname(), resid2.get_parent(), "dist =", d)
        return(d)




for res1 in residues:
    for res2 in residues:
        if res1.get_id()[1] != res2.get_id()[1]:
            aro_aro(res1, res2) # appel fonction aro_aro
