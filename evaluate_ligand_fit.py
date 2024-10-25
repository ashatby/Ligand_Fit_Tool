import fit_score_functions
import pandas as pd
import sys

if len(sys.argv)<3: 
  print("usage: python evaluate_ligand_fit.py <receptor.pdb> <ligand.pdb|docked.pdbqt>")
  sys.exit(0)

ligand_data = []

if sys.argv[2].endswith(".pdbqt"):
  # For PDBQT Files, use ligand_data_pdbqt(), parameters are:
    # proteinfile: PDB file in which the protein data is stored
    # ligandfile: PDBQT file in which the lignad data is stored
    # printprocess (Optional): Boolean variable, prints message when each ligand is being processed to monitor progress
  # Returns a 2D array of Ligand data (Name, Mean, Median, Maximum, Minimum, Q1, Q3) for each ligand

  ligand_data += (fit_score_functions.get_ligand_data_pdbqt(sys.argv[1],sys.argv[2],printprocess=True)) 


if sys.argv[2].endswith(".pdb"):
  # For PDB Files, use get_ligand_data_pdb(), parameters are:
    # proteinfile: PDB file in which the protein data is stored
    # ligandcode: HETATM Code for which ligand to process
    # hetatom_chain (Optional): If wanting to only process a specific chain of the HETATM
    # pdbname (Optional): Input for the PDB/Ligand name - By default assigns Ligand name as 4 digit code before ".pdb" in filename
    # printprocess (Optional): Boolean variable, prints message when each ligand is being processed to monitor progress
  # Returns array of Ligand data (Name, Mean, Median, Maximum, Minimum, Q1, Q3)

  ligand_data.append(fit_score_functions.get_ligand_data_pdb(sys.argv[1],sys.argv[2],hetatom_chain=sys.argv[3],pdbname="Test Ligand 1"))


# Turns 2D array of Ligand data into Pandas DataFrame 
ligand_data = fit_score_functions.dataframe(ligand_data)

# Returns only the ligands with flagged (mean > 3.8 or max > 5) values. Supports either DataFrames or List input
#print(fit_score_functions.get_flagged_values(ligand_data))

print(ligand_data.to_string())





    
