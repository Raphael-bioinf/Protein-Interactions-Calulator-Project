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
bchain = model['B']


arom_res = []
for residue in achain:
    if (residue.get_resname() == "HIS" or residue.get_resname() == "TRP" or residue.get_resname() == "TYR" or residue.get_resname() == "PHE" or residue.get_resname() == "CYS" or residue.get_resname() == "MET"):
        arom_res.append(residue)
        
        for c1 in arom_res:
            c1index = c1.get_id()[1]
            
            for c2 in arom_res:
                c2index = c2.get_id()[1]
                
                if(c1['C'] - c2['S'] < 6 and c1['C'] - c2['S'] > 5):
                    print("possible aromatic - sulfure interaction:"),
                    print (residue,c1index,"-")
                    print (residue,c2index,)
                    print (round(c1['OE2'] - c2['N'],2))
                    
                    
                    
                    



def arosulfun(resid1, resid2, dist=5.3):
    sulres = ["CYS", "MET"]
    if (resid1.get_resname() in sulres and resid2.get_resname() not in sulres):
        for atom in resid1:
            if "S" in atom.get_name():
                s_coord = atom.get_coord()
                d = dist_cal(
                    s_coord[0],
                    s_coord[1],
                    s_coord[2],
                    resid2.center_mass[0],
                    resid2.center_mass[1],
                    resid2.center_mass[2])
                if d < dist:
                    return(d)
    elif (resid2.get_resname() in sulres and resid1.get_resname() not in sulres):
        for atom in resid2:
            if "S" in atom.get_name():
                s_coord = atom.get_coord()
                d = dist_cal(
                    s_coord[0],
                    s_coord[1],
                    s_coord[2],
                    resid1.center_mass[0],
                    resid1.center_mass[1],
                    resid1.center_mass[2])
                if d < dist:
                    return(d)
