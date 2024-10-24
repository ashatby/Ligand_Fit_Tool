import KDTree
import ligand_from_pdb
import pdbqt_splicer
import pdb_1
import math
import pandas as pd


maxs = []
avgs = []

def find_fit_data(proteinfile, ligandfile):
    protein = pdb_1.read_pdb(proteinfile)
    list_of_protein_coords = list(atom.co for chain in protein for atom in chain)
    all_ligands = (pdbqt_splicer.readpdbqt(ligandfile))
    for i in range(len(all_ligands)): #each ligand in file
        ligand = all_ligands[i]
        find_ligand_fit(list_of_protein_coords,ligand,i)



def find_ligand_fit(list_of_protein_coords, ligand, itr=0):

    ligand_atom_coords = list(x.co for x in ligand)


    ribf_kdtree = KDTree.KDTree(list_of_protein_coords)


    (dsq,coord) = ribf_kdtree.nearest(ligand_atom_coords[0])
    d = math.sqrt(dsq)
    dtotal = d
    dmax = d
    dmin = d
    for i in range(1, len(ligand_atom_coords)):
        atomcords = ligand_atom_coords[i]
        (dsq,coord) = ribf_kdtree.nearest(atomcords)
        d = math.sqrt(dsq)
        # print(i,d)
        dtotal+=d
        if (d<dmin):
            dmin = d
        if(d>dmax):
            dmax = d
    davg = dtotal/len(ligand_atom_coords)
    print(f"D Average: {davg}")
    print(f'''
        LIGAND {itr}: 
        D Average: {davg}
        D Minimum: {dmin}
        D Maximum: {dmax}
    ''')

    global maxs
    global avgs
    maxs.append(dmax)
    avgs.append(davg)



def get_ligand_data_pdbqt(proteinfile, ligandfile, printprocess = False):
    # For PDBQT Files, use ligand_data_pdbqt(), parameters are:
        # proteinfile: PDB file in which the protein data is stored
        # ligandfile: PDBQT file in which the lignad data is stored
        # printprocess (Optional): Boolean variable, prints message when each ligand is being processed to monitor progress
    # Returns a 2D array of Ligand data (Name, Mean, Median, Maximum, Minimum, Q1, Q3) for each ligand

    # Example of processing ligands in RibF_vina_top100.pdbqt with protein target in ribf-af_dock.receptor_min.pdb, printing a message as each ligand is processed
        # ligand_data += (fit_score_functions.get_ligand_data_pdbqt(r"RCSB_PDBs/ribf-af_dock.receptor_min.pdb",r"RCSB_PDBs\RibF_vina_top100.pdbqt",printprocess=True)) 

    
    alldata=[]
    protein = pdb_1.read_pdb(proteinfile)
    list_of_protein_coords = list(atom.co for chain in protein for atom in chain)
    all_ligands = (pdbqt_splicer.readpdbqt(ligandfile))
    for i in range(len(all_ligands)):
        ligand = all_ligands[i]
        find_ligand_fit_pandas(list_of_protein_coords,ligand,i, alldata,pdb_name=ligand.name)
        if (printprocess):
            print("Processing Ligand",ligand.name)
    return alldata



def find_ligand_fit_pandas(list_of_protein_coords, ligand,itr, alldata, pdb_name = ''):
    data = []

    ligand_atom_coords = list(x.co for x in ligand.atoms)


    ribf_kdtree = KDTree.KDTree(list_of_protein_coords)

    fit_scores = pd.Series(dtype='float64')
    (dsq,coord) = ribf_kdtree.nearest(ligand_atom_coords[0])
    d = math.sqrt(dsq)
    fit_scores = pd.concat([fit_scores, pd.Series([d])])

    for i in range(1, len(ligand_atom_coords)):
        atomcords = ligand_atom_coords[i]
        (dsq,coord) = ribf_kdtree.nearest(atomcords)
        d = math.sqrt(dsq)
        data.append(d)
        fit_scores = pd.concat([fit_scores, pd.Series([d])])


    if (not pdb_name):
        pdb_name = "Ligand "+ str(itr+1)
    data = [pdb_name, fit_scores.mean(), fit_scores.median(), fit_scores.max(), fit_scores.min(), fit_scores.quantile(0.25), fit_scores.quantile(0.75) ]
    alldata.append(data)

    return data


def get_ligand_data_pdb(proteinfile, ligandcode,hetatom_chain='',pdbname = '',printprocess = False): 
    # For PDB Files, use get_ligand_data_pdb(), parameters are:
        # proteinfile: PDB file in which the protein data is stored
        # ligandcode: HETATM Code for which ligand to process
        # hetatom_chain (Optional): If wanting to only process a specific chain of the HETATM
        # pdbname (Optional): Input for the PDB/Ligand name - By default assigns Ligand name as 4 digit code before ".pdb" in filename
        # printprocess (Optional): Boolean variable, prints message when each ligand is being processed to monitor progress
    # Returns array of Ligand data (Name, Mean, Median, Maximum, Minimum, Q1, Q3)

    # Example of processing JRZ hetatm in 3ny8.pdb and assigning name as "Test Ligand 1"
        # ligand_data.append(fit_score_functions.get_ligand_data_pdb(r"RCSB_PDBs\3ny8.pdb","JRZ",pdbname="Test Ligand 1"))

    # Example of processing chain A of JIN hetatm in 2hzi.pdb
        # ligand_data.append(fit_score_functions.get_ligand_data_pdb(r"RCSB_PDBs\2hzi.pdb","JIN",hetatom_chain="A"))
    
    if (not pdbname):
        pdbname = proteinfile[-8:-4]
    if printprocess:
        print(f"Processing {proteinfile}")
    ligand = ligand_from_pdb.read_pdb_ligand(proteinfile,ligandcode,hetatm_chain_name=hetatom_chain,hetatms=True)
    protein = pdb_1.read_pdb(proteinfile)
    list_of_protein_coords = list(atom.co for chain in protein for atom in chain)
    ret = find_ligand_fit_pandas(list_of_protein_coords,ligand,'',[],pdb_name=pdbname)
    return ret

def dataframe(datalist):
    cols = ["Name","Mean","Median","Maximum","Minimum","Q1","Q3"]
    return pd.DataFrame(datalist,columns=cols)







