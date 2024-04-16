import bz2
import os.path
from pymatgen.core import Structure
import pandas as pd
import numpy as np
import paramiko
from crystal_analysis import Crystal
from pymatgen.analysis.bond_valence import BVAnalyzer as BVA

# SSH credentials
hostname = 'bear.dhcp.wustl.edu'
username = 'isakov'
password = 'Temp12345'
# private_key_path = '/path/to/your/private_key'

# Create a new SSH client
ssh = paramiko.SSHClient()

# Automatically add unknown hosts to the list of known hosts
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

# Load the private key for authentication
# private_key = paramiko.RSAKey.from_private_key_file(private_key_path)

def main():
    df = pd.read_csv('../data/papers/wexler/figures/data_collection_py.csv', index_col=False)

    df.replace('<null>', np.nan, inplace=True)
    df.dropna(subset=['rattle_vac_site', 'hi_sym_vac'], how='all', inplace=True)
    df.reset_index(drop=False, inplace=True)

    # # transform E_hi_sym_bulk
    df["E_hs_transformed"] = np.nan
    df["E_rb_transformed"] = np.nan

    # Conditionally calculate "E_hs_transformed"
    condition = df["E_hi_sym_bulk"].notna() & df["n_bulk_hi_sym"].notna()
    df.loc[condition, "E_hs_transformed"] = df.loc[condition, "E_hi_sym_bulk"] * df.loc[condition, "n_bulk_hi_sym"]

    # add a column to transform E_rattle_bulk to match vacancy supercell
    condition_rb = df["E_rattle_bulk"].notna() & df["n_bulk"].notna()
    df.loc[condition_rb, "E_rb_transformed"] = df.loc[condition_rb, "E_rattle_bulk"] * df.loc[condition_rb, "n_bulk"]

    # add a column for the neutral oxygen vacancy formation energy
    df["vacancy_formation_energy"] = np.nan

    # Calculate vacancy_energy based on E_hs_transformed and find minimum energy vacancy structure
    hi_sym_condition = df["E_hs_transformed"].notna()
    df.loc[hi_sym_condition, "vacancy_formation_energy"] = df["E_hi_sym_vac"] - df["E_hs_transformed"] + (-11.86096014 / 2)
    df.loc[hi_sym_condition, "group"] = df["formula"] + "_" + df["space_group"] + "_" + df["hi_sym_vac"]
    df.loc[hi_sym_condition, "poscar_name"] = df["PT_oxide"] + "_" + df["formula"] + "_" + df["space_group"] + "_" + df["hi_sym_vac"] + "_" + df["hi_sym_vac_BD"]
    df.loc[hi_sym_condition, "E_vac_min"] = df.groupby("group")["E_hi_sym_vac"].transform("min")
    df.loc[hi_sym_condition, "minimum"] = df["E_vac_min"] == df["E_hi_sym_vac"]

    # Calculate vacancy_energy based on E_rb_transformed and find mimimum energy vacancy structure
    rattle_condition = df["E_rb_transformed"].notna() & df["E_rattle_vac"].notna()
    df.loc[rattle_condition, "vacancy_formation_energy"] = df["E_rattle_vac"] - df["E_rb_transformed"] + (-11.86096014 / 2)
    df.loc[rattle_condition, "group"] = df["formula"] + "_" + df["space_group"] + "_" + df["rattle_vac_site"]
    df.loc[rattle_condition, "poscar_name"] = df["PT_oxide"] + "_" + df["formula"] + "_" + df["space_group"] + "_" + df["rattle_vac_site"] + "_" + df["rattle_vac_BD"]
    df.loc[rattle_condition, "E_vac_min"] = df.groupby("group")["E_rattle_vac"].transform("min")
    df.loc[rattle_condition, "minimum"] = df["E_rattle_vac"] == df["E_vac_min"]

    # drop any column where the vacancy structure has higher energy than the minimum
    min_E_df = df[df["minimum"]]
    # min_E_df.to_csv('looksie.csv') # three missing?

    # get poscar file to create pymatgen structure
    for _, row in min_E_df.iterrows():
        poscar_name = row["poscar_name"]
        # Connect to the remote server using SSH key authentication
        # ssh.connect(hostname=hostname, username=username, pkey=private_key)
        ssh.connect(hostname=hostname, username=username, password=password)
        file_path = os.path.join("/data3/wexler/poscars_py/", poscar_name + "_POSCAR.bz2")
        print(file_path)
        with ssh.open_sftp().open(file_path, 'rb') as f:
            poscar = bz2.decompress(f.read()).decode('utf-8')
        ssh.close()

        structure = Structure.from_str(poscar, fmt="poscar")
        try:
            structure = BVA().get_oxi_state_decorated_structure(structure)
            print('has strucuture')
            crystal = Crystal(pymatgen_structure=structure)
            print('has crystal')
        except ValueError:
            pass

if __name__ == "__main__":
    main()
