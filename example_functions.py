import fit_score_functions



ligand_data = []

# For PDBQT Files, use ligand_data_pdbqt(), parameters are:
    # proteinfile: PDB file in which the protein data is stored
    # ligandfile: PDBQT file in which the lignad data is stored
    # printprocess (Optional): Boolean variable, prints message when each ligand is being processed to monitor progress
# Returns a 2D array of Ligand data (Name, Mean, Median, Maximum, Minimum, Q1, Q3) for each ligand

# Example of processing ligands in RibF_vina_top100.pdbqt with protein target in ribf-af_dock.receptor_min.pdb, printing a message as each ligand is processed
ligand_data += (fit_score_functions.get_ligand_data_pdbqt(r"RCSB_PDBs/ribf-af_dock.receptor_min.pdb",r"RCSB_PDBs\RibF_vina_top100.pdbqt",printprocess=True)) 




# For PDB Files, use get_ligand_data_pdb(), parameters are:
    # proteinfile: PDB file in which the protein data is stored
    # ligandcode: HETATM Code for which ligand to process
    # hetatom_chain (Optional): If wanting to only process a specific chain of the HETATM
    # pdbname (Optional): Input for the PDB/Ligand name - By default assigns Ligand name as 4 digit code before ".pdb" in filename
    # printprocess (Optional): Boolean variable, prints message when each ligand is being processed to monitor progress
# Returns array of Ligand data (Name, Mean, Median, Maximum, Minimum, Q1, Q3)

# Example of processing JRZ hetatm in 3ny8.pdb and assigning name as "Test Ligand 1"
ligand_data.append(fit_score_functions.get_ligand_data_pdb(r"RCSB_PDBs\3ny8.pdb","JRZ",pdbname="Test Ligand 1"))

# Example of processing chain A of JIN hetatm in 2hzi.pdb
ligand_data.append(fit_score_functions.get_ligand_data_pdb(r"RCSB_PDBs\2hzi.pdb","JIN",hetatom_chain="A"))

# Turns 2D array of Ligand data into Pandas DataFrame 
ligand_data = fit_score_functions.dataframe(ligand_data)

print(ligand_data)





    