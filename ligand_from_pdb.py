import sys
import math

# returns a list of chains, which are a list of ATOM records
import KDTree
import pdbqt_splicer

def read_pdb(filename,hetatms=False,hydrogens=False):
  file = open(filename,"r")
  lines = file.readlines()
  file.close()
  return parse_pdb(lines,hetatms)

def read_pdb_ligand(filename,hetatm_name,hetatm_chain_name = '',hetatms=False,hydrogens=False):
  file = open(filename,"r")
  lines = file.readlines()
  file.close()
  protein = parse_pdb(lines,hetatms)

#   for chain in protein:
#     for atom in chain.atoms:
#       if atom.is_hetatm:
#         print(atom.chn)


#   hetatms = [chain for chain in protein if all(atom.is_hetatm for atom in chain.atoms)]
#   ligand = [atom for chain in hetatms for atom in chain.atoms if atom.rnam == hetatm_name]

  hetatms = [atom for chain in protein for atom in chain.atoms if atom.is_hetatm]
  if (hetatm_chain_name):
    ligand = [atom for atom in hetatms if ((atom.rnam == hetatm_name) and (atom.chn == hetatm_chain_name))]
  else:
    ligand = [atom for atom in hetatms if atom.rnam == hetatm_name]
    # print(f"ligand het chain is '{ligand[0].chn}' for filename {filename}")

    
  return ligand

 



def parse_pdb(lines,hetatms=False,hydrogens=False):
  chains = []
  atoms = []
  for line in lines:
    if line.find("ATOM")==0:
      if line[13]!='H' or hydrogens==True:
        if line[16]==' ' or line[16]=='A': # only first atom if multiple altLoc
          atoms.append(atom_rec(line))
    if line.find("HETATM")==0 and hetatms==True: atoms.append(atom_rec(line, True))
    if line.find("TER")==0 or line.find("ENDMDL")==0: 
      chains.append(protein_chain(atoms))
      atoms = []
  if len(atoms)>0: chains.append(protein_chain(atoms))
  return chains #a 2d array of atoms

# parse ATOM records into (anum,anam,rnum,rnam,chn,x,y,z,occ,bf)
class protein_chain:
  def __init__(self,iatoms):
    self.atoms = iatoms
    # self.name = iatoms[0].rnam
    self.hetmolec = iatoms[0].is_hetatm
    # self.hetmolec = ihetatm

    def return_hetatm_chain(self, hetatm_name):
      return list(atom for atom in self.atoms if ((atom.rnam == hetatm_name) and (atom.is_hetatm))) 



class atom_rec:
  def __init__(self,line,hetatm = False):
    self.line = line
    self.anum = int(line[6:11])
    self.anam = line[12:16]
    self.rnam = line[17:20].strip()
    self.chn = line[21:22]
    self.rnum = int(line[22:26])
    self.x = float(line[30:38])
    self.y = float(line[38:46])
    self.z = float(line[46:54])
    self.co = [self.x,self.y,self.z]
    self.occ = float(line[54:60])
    self.bf = float(line[60:66])
    self.is_hetatm = (line[0:6] == "HETATM")


  def toString(self):
    return "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" % \
      (self.anum,self.anam,self.rnam,self.chn,self.rnum%10000, \
       self.co[0],self.co[1],self.co[2],self.occ,self.bf)
  def copy(self): return atom_rec(self.toString())


# # test = sys.argv[1] # "4HHB.pdb"
# test = "RCSB_PDBs/4ph9.pdb"
# pdb = read_pdb(test)
# first_atom_chain = pdb[0]
# first_atom = first_atom_chain[0]
# next_atom = first_atom_chain[2]
# other_atoms = first_atom_chain[5:1:-1]
# # print ("num_chains("+test+")="+str((pdb)[0][0].line))
# print(f'''
#   self.line = {first_atom.line}
#   self.anum = {first_atom.anum}
#   self.anam = {first_atom.anam}
#   self.rnam = {first_atom.rnam}
#   self.chn = {first_atom.chn}
#   self.rnum = {first_atom.rnum}
#   self.x = {first_atom.x}
#   self.y = {first_atom.y}
#   self.z = {first_atom.z}
#   self.co = {first_atom.co}
#   self.occ = {first_atom.occ}
#   self.bf = {first_atom.bf}   

# -------

#   self.line = {next_atom.line}
#   self.anum = {next_atom.anum}
#   self.anam = {next_atom.anam}
#   self.rnam = {next_atom.rnam}
#   self.chn = {next_atom.chn}
#   self.rnum = {next_atom.rnum}
#   self.x = {next_atom.x}
#   self.y = {next_atom.y}
#   self.z = {next_atom.z}
#   self.co = {next_atom.co}
#   self.occ = {next_atom.occ}
#   self.bf = {next_atom.bf} 

# -------

# THE TWO COORDS ARE: {first_atom.co} and {next_atom.co}
#       ''')


# print(f"COORDS FOR ORIGINAL ATOM: {first_atom.co}")

# other_atom_coords = []
# for i in range (len(other_atoms)):
#   atom = other_atoms[i]
#   print(f"COORDS FOR ATOM {i}: {atom.co}")
#   other_atom_coords.append(atom.co)


# # example of how to use this:
# #  kdtree = KDTree([x.co for x in receptor])
# #  (dsq,coord) = kdtree.nearest(prot[j].co)
# #  j = kdtree.index_of(coord)


# list_of_atoms = [first_atom.co, next_atom.co]
# print(list_of_atoms)
# print(other_atom_coords)

# kdtree = KDTree.KDTree(other_atom_coords)
# (dsq,coord) = kdtree.nearest(first_atom.co)
# j = kdtree.index_of(coord)
# d = math.sqrt(dsq)
# print(j,d)


