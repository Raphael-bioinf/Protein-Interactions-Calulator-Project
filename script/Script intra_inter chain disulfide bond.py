prot_id = "1eej.pdb"
prot_file = "1eej.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]
achain = model['A']
bchain = model['B']
    
cysresiduesA = []
cysresiduesB = []
for residue in achain:
    if residue.get_resname() == "CYS":
        cysresiduesA.append(residue)
for c1 in cysresiduesA:
    c1index = c1.get_id()[1]
    for c2 in cysresiduesA:
        c2index = c2.get_id()[1]
        if(c1index != c2index):
            if(c1['SG'] - c2['SG'])< 2.2:
                print("possible di-sulfide bond in chain A:"),
                print ("Cys",c1index,"-")
                print ("Cys",c2index,)
                print (round(c1['SG'] - c2['SG'],2))
            
for residue in bchain:
    if residue.get_resname() == "CYS":
        cysresiduesB.append(residue)
for c3 in cysresiduesB:
    c3index = c3.get_id()[1]
    for c4 in cysresiduesB:
        c4index = c4.get_id()[1]
        if(c3index != c4index):
            if (c3['SG'] - c4['SG'])< 2.2:
                print("possible di-sulfide bond in chain B:"),
                print ("Cys",c3index,"-")
                print ("Cys",c4index,)
                print (round(c3['SG'] - c4['SG'],2))
                
for c1 in cysresiduesA:
    c1index = c1.get_id()[1]
    for c2 in cysresiduesA:
        c2index = c2.get_id()[1]
        for c3 in cysresiduesB:
                c3index = c3.get_id()[1]
                for c4 in cysresiduesB:
                        c4index = c4.get_id()[1]
                        if(c1index != c2index):
                            if(c1['SG'] - c2['SG'])< 2.2:
                                print("possible di-sulfide bond in chain A:"),
                                print ("Cys",c1index,"-")
                                print ("Cys",c2index,)
                                print (round(c1['SG'] - c2['SG'],2))
                        if(c3index != c4index):
                            if (c3['SG'] - c4['SG'])< 2.2:
                                print("possible di-sulfide bond in chain B:"),
                                print ("Cys",c3index,"-")
                                print ("Cys",c4index,)
                                print (round(c3['SG'] - c4['SG'],2))
                        if(c1['SG'] - c3['SG'])< 2.2:
                            print("possible di-sulfide bond between chain A and B:"),
                            print ("Cys",c1index,"-")
                            print ("Cys",c3index,)
                            print (round(c1['SG'] - c3['SG'],2))
                        if(c2['SG'] - c4['SG'])< 2.2:
                            print("possible di-sulfide bond between chain A and B:"),
                            print ("Cys",c2index,"-")
                            print ("Cys",c4index,)
                            print (round(c2['SG'] - c4['SG'],2))
                        if(c1['SG'] - c4['SG'])< 2.2:
                            print("possible di-sulfide bond between chain A and B:"),
                            print ("Cys",c1index,"-")
                            print ("Cys",c4index,)
                            print (round(c1['SG'] - c4['SG'],2))
                        if(c2['SG'] - c3['SG'])< 2.2:
                            print("possible di-sulfide bond between chain A and B:"),
                            print ("Cys",c2index,"-")
                            print ("Cys",c3index,)
                            print (round(c2['SG'] - c3['SG'],2))






