# Ligand_Fit_Tool
 This is a tool made by Anthony Shatby under Dr. Thomas Ioerger to determine ligand fit into a protein.
 Please see (or run) **example_functions.py** to see functionality and syntax.

Notes:
 - if the pdbqt file contains multiple ligands, it is assumed they are docked by Autodock-VINA, which uses 'REMARK  NAME' to indicate compound id's
 - hydrogen atoms on the protein and ligand are ignored
 - For a given PDB file, ensure that the ligand in interest has all nearby protein residues present
   - For example, in PDB 3L5D, when looking at ligand BDV, chain B was used exclusively
   - This was because chain A of the protein (next to chain A of the BDV Ligand) was missing residues 131-136
   - The true closest protein atom was found on TYR 132, which was missing from chain A
   - This created an incorrect, artificially high distance score for this atom, and subsequently the entire ligand
   - This was avoided by looking exclusively analyzing chain B, which had all protein residues properly crystalized

Here is the usage:

```
  # python evaluate_ligand_fit.py <receptor.pdb> <ligands.pdbqt>
```

It prints out a table with various metrics on each ligand-protein interface.
It starts by finding the closest protein atom for each ligand atom.
The metrics include the min, mean, max of these distances, along with 
Q1 and Q3 (25% and 75% quartiles).

If a PDB file of a complex is provided (containing the ligand), it can be indicated by giving the 
code ('residue name') of the ligand and chain id.  It is assumed the ligand coordinates
are given in HETATM records.  For example
 
```
  # python evaluate_ligand_fit.py 8J8J.pdb PRP D
```
