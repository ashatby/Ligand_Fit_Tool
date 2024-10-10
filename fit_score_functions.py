import KDTree
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
    for i in range(len(all_ligands)): #each ligand in top 100 file
        ligand = all_ligands[i]
        find_ligand_fit(list_of_protein_coords,ligand,i)
        print("Processing Ligand",i)



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



def find_fit_data_pandas(proteinfile, ligandfile):
    alldata=[]
    protein = pdb_1.read_pdb(proteinfile)
    list_of_protein_coords = list(atom.co for chain in protein for atom in chain)
    all_ligands = (pdbqt_splicer.readpdbqt(ligandfile))
    for i in range(len(all_ligands)):
        ligand = all_ligands[i]
        find_ligand_fit_pandas(list_of_protein_coords,ligand,i, alldata)
        print("Processing Ligand",i)
    df = pd.DataFrame(alldata)
    df.columns = ["Name","Mean","Median","Max","Min","Q1","Q3"]
    print(df)
    df.to_csv("liganddata_vina.csv")


def find_ligand_fit_pandas(list_of_protein_coords, ligand,itr, alldata, pdb_name = ''):
    data = []

    ligand_atom_coords = list(x.co for x in ligand)


    ribf_kdtree = KDTree.KDTree(list_of_protein_coords)

    fit_scores = pd.Series(dtype='float64')
    (dsq,coord) = ribf_kdtree.nearest(ligand_atom_coords[0])
    d = math.sqrt(dsq)
    # data.append(d)
    fit_scores = pd.concat([fit_scores, pd.Series([d])])

    for i in range(1, len(ligand_atom_coords)):
        atomcords = ligand_atom_coords[i]
        (dsq,coord) = ribf_kdtree.nearest(atomcords)
        d = math.sqrt(dsq)
        data.append(d)
        fit_scores = pd.concat([fit_scores, pd.Series([d])])
#     print(f'''  
#                 Mean: {fit_scores.mean()}
#                 Median: {fit_scores.median()}
#                 Maximum: {fit_scores.max()}
#                 Minimum: {fit_scores.min()}
#                 Q1: {fit_scores.quantile(0.25)}
#                 Q3: {fit_scores.quantile(0.75)}
# ''')
    data = [pdb_name, fit_scores.mean(), fit_scores.median(), fit_scores.max(), fit_scores.min(), fit_scores.quantile(0.25), fit_scores.quantile(0.75) ]
    # print(data)
    return data
    # alldata.append(data)


proteinfile = "RCSB_PDBs/ribf-af_dock.receptor_min.pdb"
ligandsfile = "RCSB_PDBs/RibF_vina_top100.pdbqt"
# ligandsfile = r"C:\Users\shatb\Downloads\Research\RCSB_PDBs\RibF_SACC_hits103.pdbqt"

# find_fit_data_pandas(proteinfile, ligandsfile)
# find_fit_data_pandas(proteinfile, ligandsfile)

# print("Maximum average distance:",max(avgs))
# print("Minimum average distance:",min(avgs))




