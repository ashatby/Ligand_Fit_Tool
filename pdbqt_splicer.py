import re

class ligand_atom:
  def __init__(self,line):
    atom_attributes = ((lambda x: [i.strip() for i in x.split()])(line))
    self.line = line
    self.anum = int(atom_attributes[1])
    self.anam = atom_attributes[2]
    self.rnam = atom_attributes[3]
    self.rnum = int(atom_attributes[4])
    self.x = float(atom_attributes[5])
    self.y = float(atom_attributes[6])
    self.z = float(atom_attributes[7])
    self.co = [self.x,self.y,self.z]
    self.occ = float(atom_attributes[8])
    self.bf = float(atom_attributes[9])
    self.partialcharge = float(atom_attributes[10])
    self.element = atom_attributes[11]
#   def toString(self):
#     return "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" % \
#       (self.anum,self.anam,self.rnam,self.chn,self.rnum%10000, \
#        self.co[0],self.co[1],self.co[2],self.occ,self.bf)
#   def copy(self): return atom_rec(self.toString())

def readpdbqt(filename):
    file = open(filename,"r")
    lines = file.readlines()
    file.close()
    return parse_pdbqt(lines)

def parse_pdbqt(lines): #return a 2d array of ligands, atoms in ligand
    ligands = []
    atoms = []
    
    molecules = splice_into_molecules_from_textlines(lines)
    # print(len(molecules), "is molecs")
    # print(molecules[0])
    # print("--------------")
    # print(molecules[1])


    for molecule in molecules:
        for line in molecule.splitlines():
            if line.startswith("ATOM"):
                newatom = extract_atom_attributes_from_line(line)
                # print("appended newatom")
                atoms.append(newatom)
        ligands.append(atoms)
        # print("appended ligands is")
        # print(atoms)
        atoms = []
    # print("ligands")

    # print(ligands[0][0].co)
    # print(ligands[1][0].co)

    return ligands

    # for line in lines:
    #     if line.find("ATOM")==0:
    #         if line[13]!='H' or hydrogens==True:
    #             if line[16]==' ' or line[16]=='A': # only first atom if multiple altLoc
    #                 atoms.append(atom_rec(line))
    #     if line.find("HETATM")==0 and hetatms==True: atoms.append(atom_rec(line))
    #     if line.find("TER")==0 or line.find("ENDMDL")==0: 
    #             chains.append(atoms)
    #             atoms = []
    # if len(atoms)>0: chains.append(atoms)
    # return chains #a 2d array of atoms

def extract_atom_attributes_from_line(line):
    if line.startswith("ATOM"):
        newatom = ligand_atom(line)
        # print(newatom)
        return newatom
    else:
        return "ERROR IN INPUT LINE"


def splice_into_molecules_from_textlines(text): #returns list of raw lines from pdbqt file, split by molecule
    # print(text)
    stringtext = "".join(text)
    start_word = "MODEL"
    end_word = "ENDMDL"

    # Regular expression to find the content between start_word and end_word
    pattern = re.compile(f'{re.escape(start_word)}(.*?){re.escape(end_word)}', re.DOTALL)
    substrings = pattern.findall(stringtext)
    return substrings

# start_word = "MODEL"
# end_word = "ENDMDL"

# # Regular expression to find the content between start_word and end_word
# pattern = re.compile(f'{re.escape(start_word)}(.*?){re.escape(end_word)}', re.DOTALL)
# substrings = pattern.findall(text)

# substrings = get_pdbqt_molecules_from_file("RCSB_PDBs/RibF_vina_top100.pdbqt")

# print((substrings))
# print(len(substrings))

# for i in range(len(substrings)):
#     substring = substrings[i]
#     print(f'''
# MOLECULE: {i+1}
#        {substring}             
# END OF MOLECULE: {i+1}

# ''')
    
# substring = substrings[0]
# for line in substring.splitlines():
#     if line.startswith("ATOM"):
#         print(line)
#     print("DOESNT HAVE ATOM")




# readpdbqt("RCSB_PDBs/RibF_SACC_hits103.pdbqt")
# readpdbqt("RCSB_PDBs/RibF_vina_top100.pdbqt")

    

