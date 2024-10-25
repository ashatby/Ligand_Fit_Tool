# Ligand_Fit_Tool
 This is a tool made by Anthony Shatby under Dr. Thomas Ioerger to determine ligand fit into a protein.
 Please see (or run) **example_functions.py** to see functionality and syntax.

Notes:
 - if the pdbqt file contains multiple ligands, it is assumed they are docked by Autodock-VINA, which uses 'REMARK  NAME' to indicate compound id's
 - hydrogren atoms on the protein and ligand are ignored

Here is the usage:

```
  # python evaluate_ligand_fit.py <receptor.pdb> <ligands.pdbqt>
```

It prints out a table with various metrics on the ligand-protein interface.
It starts by finding the closest protein atom for each ligand atom.
The metrics include the min, mean, max of these distances, along with 
Q1 and Q3 (25% and 75% quartiles).

If a PDB file of a complex is provided (containing the ligand), it can be indicated by giving the 
code ('residue name') of the ligand and chain id.  It is assumed the ligand coordinates
are given in HETATM records.
 