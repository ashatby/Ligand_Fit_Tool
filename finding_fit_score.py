import KDTree
import pdbqt_splicer
import pdb_1
import math



proteinfile = "RCSB_PDBs/ribf-af_dock.receptor_min.pdb"

#########################################################
# ribf = pdb_1.read_pdb(proteinfile)
# list_of_atom_coords = (list(atom.co for atom in ribf[0]))
# print("lofff",len(list_of_atom_coords))
#########################################################


#########################################################
ribf = pdb_1.read_pdb(proteinfile)
list_of_atom_coords = list(atom.co for chain in ribf for atom in chain)
#########################################################



all_ligands = (pdbqt_splicer.readpdbqt("RCSB_PDBs/RibF_vina_top100.pdbqt"))
first_ligand = all_ligands[0]
first_ligand_atom = all_ligands[0][0]
first_ligand_atom_coords = list(x.co for x in first_ligand)

# print((list(ribf[0])))
prots = ((list(ribf[0])))

# print(type(list(ribf)[0][0]))

print(type(first_ligand[0]))

# for i in range(len(first_ligand_atom_coords)):
#     print(first_ligand_atom_coords[i])


ribf_kdtree = KDTree.KDTree(list_of_atom_coords)
(dsq,coord) = ribf_kdtree.nearest(first_ligand_atom.co)
j = ribf_kdtree.index_of(coord)
d = math.sqrt(dsq)
print(j,d)




ribf_kdtree = KDTree.KDTree(list_of_atom_coords)


(dsq,coord) = ribf_kdtree.nearest(first_ligand_atom.co)
d = math.sqrt(dsq)
dtotal = d
dmax = d
dmin = d
for i in range(1, len(first_ligand_atom_coords)):
    atomcords = first_ligand_atom_coords[i]
    (dsq,coord) = ribf_kdtree.nearest(atomcords)
    d = math.sqrt(dsq)
    print(i,d)
    dtotal+=d
    if (d<dmin):
        dmin = d
    if(d>dmax):
        dmax = d
davg = dtotal/len(first_ligand_atom_coords)

print(f"D Average: {davg}")
print(f'''
    D Average: {davg}
    D Minimum: {dmin}
    D Maximum: {dmax}
''')


