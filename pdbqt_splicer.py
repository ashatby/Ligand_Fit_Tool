import re

class ligand_molec:
  def __init__(self,lines):
    self.atoms = []
    self.name = ''
    for line in lines:
            # Look for the line that starts with 'REMARK file:'
            if line.startswith("REMARK file:"):
                # Extract the name part (after "REMARK file:") and strip any trailing spaces/newlines
                name_with_extension = line.split("REMARK file:")[1].strip()
                # Remove the '.pdbqt' extension
                if name_with_extension.endswith(".pdbqt"):
                    self.name = name_with_extension[:-6]  # Remove last 6 characters (".pdbqt")
                break



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

    return ligands



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



    

