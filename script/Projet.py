#!/usr/bin/env python
# coding: utf-8
# In[]:

# Programme permettant d'accéder aux données et aux coordonnées du fichier PDB


prot_id = "5AGY.pdb"
prot_file = "../data/5AGY.pdb"



amino = {
    "hydrophobic": [
        "ALA",
        "VAL",
        "LEU",
        "ILE",
        "MET",
        "PHE",
        "TRP",
        "PRO",
        "TYR"],
    "disulphide": ["CYS"],
    "ionic": [
        "ARG",
        "LYS",
        "HIS",
        "ASP",
        "GLU"],
    "aroaro": [
        "PHE",
        "TRP",
        "TYR"],
    "arosul": [
        "PHE",
        "TRP",
        "TYR",
        "CYS",
        "MET"],
    "cationpi": [
        "LYS",
        "ARG",
        "PHE",
        "TRP",
        "TYR"],
    "all": [
        "ALA",
        "GLY",
        "ILE",
        "LEU",
        "PRO",
        "VAL",
        "PHE",
        "TRP",
        "TYR",
        "ASP",
        "GLU",
        "ARG",
        "HIS",
        "LYS",
        "SER",
        "THR",
        "CYS",
        "MET",
        "ASN",
        "GLN"]}

from math import sqrt
from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]


for chain in model:
    for residue in chain:
        for atom in residue:
            print (chain, residue, atom, atom.get_coord())







