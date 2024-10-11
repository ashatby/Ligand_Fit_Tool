import ligand_from_pdb
import fit_score_functions
import pdb_1
import pandas as pd

def get_ligand_data(proteinfile, ligandcode,alldata,hetatom_chain=''):
    print(f"Processing {proteinfile}")
    ligand = ligand_from_pdb.read_pdb_ligand(proteinfile,ligandcode,hetatm_chain_name=hetatom_chain,hetatms=True)
    protein = pdb_1.read_pdb(proteinfile)
    list_of_protein_coords = list(atom.co for chain in protein for atom in chain)
    # print(fit_score_functions.find_ligand_fit_pandas(list_of_protein_coords,ligand,'',[],pdb_name=proteinfile[-8:-4]))
    alldata.append(fit_score_functions.find_ligand_fit_pandas(list_of_protein_coords,ligand,'',[],pdb_name=proteinfile[-8:-4]))



def make_pdb_file(pdb):
    return f"RCSB_PDBs/{pdb}.pdb"




pdbs_to_code = {"3eml": "ZMA",
                "3bkl": "KAW",
                "1e66": "HUX",
                "2e1w": "FR6",
                "2oi0": "283",
                "3cqw": "CQW",
                "2am9": "TES",
                "1bcd": "FMS",
                "2cnk": "MY2",
                "3ny8": "JRZ",
                "2hv5": "ZST",
                "1h00": "FCP"
                }

pdbs_to_code_chain_a = {
                "2hzi": "JIN",
                "2vt4": "P32",
                "3d0e": "G93",
                "1s3b": "RMA",
                "3l5d": "BDV",
                "3d4q": "SM5"
                }


# ligand = ligand_from_pdb.read_pdb_ligand(proteinfile,"ZMA",hetatms=True)
# protein = pdb_1.read_pdb(proteinfile)
# list_of_protein_coords = list(atom.co for chain in protein for atom in chain)

# # fit_score_functions.find_ligand_fit_pandas(list_of_protein_coords,ligand,0,[],pdb_name=proteinfile[-8:-4])
cols = (["PDB File","Mean","Median","Max","Min","Q1","Q3"])
# cols = (["PDB File","Ligand","Chain","Mean","Median","Max","Min","Q1","Q3"])
# print_ligand_data(_3emlfile, "ZMA")
# print_ligand_data(_3bklfile, "KAW")
alldata=[]
for pdb,code in pdbs_to_code.items():
    get_ligand_data(make_pdb_file(pdb), code,alldata)

# df = pd.DataFrame(alldata, columns=cols)

# print(df)
for pdb,code in pdbs_to_code_chain_a.items():
    get_ligand_data(make_pdb_file(pdb), code,alldata,"A")

# # pdb = "2hzi"
# (get_ligand_data(make_pdb_file("2hzi"), "JIN",alldata,"A"))
# (get_ligand_data(make_pdb_file("2vt4"), "P32",alldata,"A"))
# (get_ligand_data(make_pdb_file("3d0e"), "G93",alldata,"A"))
# (get_ligand_data(make_pdb_file("1s3b"), "RMA",alldata,"A"))
# (get_ligand_data(make_pdb_file("3l5d"), "BDV",alldata,"A"))
# (get_ligand_data(make_pdb_file("3d4q"), "SM5",alldata,"A"))

print(alldata)

df = pd.DataFrame(alldata, columns=cols)

print(df)

df.to_csv("ligand_top_20_dude.csv")



'''
"2hzi" - JIN A
"2vt4" - P32 A
"3d0e" - G93 A
"1s3b" - RMA A
"3l5d" - BDV A
"3d4q" - SM5 A

'''



'''
3eml ZMA
3bkl KAW
1e66 HUX
2e1w FR6
2oi0 283
3cqw CQW
2am9 TES
1bcd FMS
2cnk MY2?
1h00 FCP
2hv5 ZST
1h00 FCP

'''



