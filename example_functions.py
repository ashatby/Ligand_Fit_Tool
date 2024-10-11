import pdb_1
import ligand_from_pdb
import fit_score_functions
import pandas as pd
import KDTree
import pdbqt_splicer
import math

# print(pdb_1.read_pdb(r"RCSB_PDBs/1e66.pdb")) #returns a list of chains. each chain is a list of atoms. 2D array
alldata=[]
# alldata.append(fit_score_functions.get_ligand_data(r"RCSB_PDBs\3ny8.pdb","JRZ",alldata))
# alldata.append(fit_score_functions.get_ligand_data(r"RCSB_PDBs\2hzi.pdb","JIN",hetatom_chain="A"))
for row in alldata:
    print(row)


alldata=[]



alldata = (fit_score_functions.find_fit_data_pandas(r"RCSB_PDBs/ribf-af_dock.receptor_min.pdb",r"RCSB_PDBs\RibF_vina_top100.pdbqt")) #put in protein (PDB) and ligand (PDBQT) files. returns 2d array of ligand data 

for row in alldata:
    print(row)