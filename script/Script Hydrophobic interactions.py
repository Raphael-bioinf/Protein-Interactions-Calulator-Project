
#PHYDROPHOBIC INTERACTIONS:The following residues are considered to participate in interactions if they fall within 5Å range.
#ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, TYR.

prot_id = "1eej.pdb"
prot_file = "1eej.pdb"

prot_id = "1gai.pdb"
prot_file = "1gai.pdb"

from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure(prot_id, prot_file)
model=structure[0]

residuesA = []
#On selectionne les acides aminés impliqués dans les intercations hydrophobes
for chain in model:
    hydrophobicAA=["ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","TYR"]
    for residue in chain:
        if residue.get_resname() in hydrophobicAA:
            residuesA.append(residue)
#Les atomes présents dans la chaine principale des acides aminés sont considérés comme ne participant pas à ces interactions
for r1 in residuesA:
    cresidues =['CA''O','N','H','C'] 
    for atom in r1:
            r1name=atom.get_name()
            if r1name not in cresidues:
                #On récupère le numéro affecté à chaque residue (acide aminé)
                    r1index = r1.get_id()[1]
    for r2 in residuesA:
        for atom in r2:
            r2name=atom.get_name()
            if r2name not in cresidues:
                r2index = r2.get_id()[1]
                #On effectue la comparaison pour la  distance entre deux atomes
                if(r1index != r2index):
                    try:
                        distance = r1[r1name] - r2[r2name]
                    except KeyError:
                        continue
                    if(r1[r1name] - r2[r2name])< 5:
                        print("possible hydrophobic interactions in chain:",residue.get_parent()),
                        print (r1.get_resname(),r1index,"-",r2.get_resname(),r2index,)
                        print (round(r1[r1name] - r2[r2name],2))
