import bz2
import os.path

from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core import Structure, Composition
import pandas as pd
import numpy as np
import paramiko
from crystal_analysis import Crystal
from pymatgen.analysis.bond_valence import BVAnalyzer as BVA
from pymatgen.transformations.standard_transformations import OxidationStateDecorationTransformation

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
    df = pd.read_csv('../data/papers/wexler/figures/data_collection_py0.csv', index_col=False)

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
    hi_sym_condition = df["E_hs_transformed"].notna() & df["E_hi_sym_vac"].notna()
    df.loc[hi_sym_condition, "vacancy_formation_energy"] = df["E_hi_sym_vac"] - df["E_hs_transformed"] + (-11.86096014 / 2)
    df.loc[hi_sym_condition, "Eg"] = df["bandgap_hi_sym_vac"]
    df.loc[hi_sym_condition, "bulk"] = "higher_symmetry"
    df.loc[hi_sym_condition, "E_bulk"] = df["E_hs_transformed"]
    df.loc[hi_sym_condition, "site"] = df["hi_sym_vac"]
    df.loc[hi_sym_condition, "group"] = df["formula"] + "_" + df["space_group"] + "_" + df["hi_sym_vac"]
    df.loc[hi_sym_condition, "poscar_name"] = df["PT_oxide"] + "_" + df["formula"] + "_" + df["space_group"] + "_" + df["bulk"]
    # df.loc[hi_sym_condition, "vacancy_poscar_name"] = df["PT_oxide"] + "_" + df["formula"] + "_" + df["space_group"] + "_" + df["hi_sym_vac"] + "_" + df["hi_sym_vac_BD"]
    df.loc[hi_sym_condition, "E_vac_min"] = df.groupby("group")["E_hi_sym_vac"].transform("min")
    df.loc[hi_sym_condition, "minimum"] = df["E_vac_min"] == df["E_hi_sym_vac"]

    # Calculate vacancy_energy based on E_rb_transformed and find mimimum energy vacancy structure
    rattle_condition = df["E_rb_transformed"].notna() & df["E_rattle_vac"].notna()
    df.loc[rattle_condition, "vacancy_formation_energy"] = df["E_rattle_vac"] - df["E_rb_transformed"] + (-11.86096014 / 2)
    df.loc[rattle_condition, "Eg"] = df["bandgap_rattle_vac"]
    df.loc[rattle_condition, "bulk"] = "rattle"
    df.loc[rattle_condition, "E_bulk"] = df["E_rb_transformed"]
    df.loc[rattle_condition, "site"] = df["rattle_vac_site"]
    df.loc[rattle_condition, "group"] = df["formula"] + "_" + df["space_group"] + "_" + df["rattle_vac_site"]
    df.loc[rattle_condition, "poscar_name"] = df["PT_oxide"] + "_" + df["formula"] + "_" + df["space_group"]
    # df.loc[rattle_condition, "vacancy_poscar_name"] = df["PT_oxide"] + "_" + df["formula"] + "_" + df["space_group"] + "_" + df["rattle_vac_site"] + "_" + df["rattle_vac_BD"]
    df.loc[rattle_condition, "E_vac_min"] = df.groupby("group")["E_rattle_vac"].transform("min")
    df.loc[rattle_condition, "minimum"] = df["E_rattle_vac"] == df["E_vac_min"]

    rerun_condition = df["E_rb_transformed"].notna() & df["E_rattle_vac"].notna() & df["vac_rerun"].notna()
    df.loc[rerun_condition, "vacancy_formation_energy"] = df["E_vac_rerun"] - df["E_rb_transformed"] + (-11.86096014 / 2)
    hs_rerun_condition = df["E_hs_transformed"].notna() & df["E_hi_sym_vac"].notna() & df["vac_rerun"].notna()
    df.loc[hs_rerun_condition, "vacancy_formation_energy"] = df["E_vac_rerun"] - df["E_hs_transformed"] + (-11.86096014 / 2)

    # drop any column where the vacancy structure has higher energy than the minimum
    min_E_df = df[df["minimum"]]
    # min_E_df.to_csv('code_minE_agrees.csv')

    # keep only the calculations rerun at a denser grid (note these match exactly with those left from the minimums)
    min_E_df = df[df["vac_rerun"].notna()]
    # min_E_df.to_csv('rerun_vac_Ev.csv')

    # remove rattle vacancies that exist when rattle bulk is at higher energy than hi_sym
    # min_E_df["low_E_structure"] = min_E_df.groupby("group")["E_vac_min"].transform("min")
    # min_E_df["stable_structure"] = min_E_df["E_vac_min"] == min_E_df["low_E_structure"]
    min_E_df["low_E_structure"] = min_E_df.groupby(["formula", "space_group"])["E_bulk"].transform("min")
    min_E_df["stable_structure"] = min_E_df["E_bulk"] == min_E_df["low_E_structure"]

    min_E_df = min_E_df[min_E_df["stable_structure"]]
    # min_E_df.to_csv("wexler_full_processed.csv"
    df_processed = min_E_df[["PT_oxide", "formula", "space_group", "MP_ehull", "MP_Eg", "E_vac_rerun",
                            "vacancy_formation_energy", "Eg", "E_bulk", "group", "poscar_name", "E_vac_min"]].reset_index(drop=True)
    df_processed.to_csv("wexler_processed.csv")

    df_cf = pd.DataFrame()
    # get poscar file to create pymatgen structure
    for _, row in min_E_df.iterrows():
        poscar_name = row["poscar_name"]
        print(row["group"])
        formula = row["formula"]
        space_group = row["space_group"]
        # Connect to the remote server using SSH key authentication
        # ssh.connect(hostname=hostname, username=username, pkey=private_key)
        ssh.connect(hostname=hostname, username=username, password=password)
        file_path = os.path.join("/data3/wexler/poscars_py/", poscar_name + "_POSCAR.bz2")
        with ssh.open_sftp().open(file_path, 'rb') as f:
            poscar = bz2.decompress(f.read()).decode('utf-8')

        ssh.close()
        # crystal = Crystal(poscar_string=poscar)
        # print(crystal)

        structure = Structure.from_str(poscar, fmt="poscar")
        index = int(row["site"].split("s")[1].replace("_0", ""))
        # print(index)

        composition = Composition(formula)
        oxi_dict = Composition.oxi_state_guesses(composition, target_charge=0)[0]
        oxi_states = OxidationStateDecorationTransformation(oxi_dict)
        structure = oxi_states.apply_transformation(structure)
        # print(structure)
        # print(structure[index].frac_coords)
        try:
            crystal = Crystal(pymatgen_structure=structure, n=index)
            # crystal = Crystal(pymatgen_structure=structure, n=index,
            #                   nn_finder=CrystalNN(weighted_cn=True, cation_anion=True), use_weights=True)

            CN = crystal.cn_dicts
            Eb = crystal.bond_dissociation_enthalpies
            Vr = crystal.reduction_potentials

            Eb_sum = []
            for CN_dict, Eb_dict in zip(CN, Eb):
                CN_array = np.array(list(CN_dict.values()))
                Eb_array = np.array(list(Eb_dict.values()))
                Eb_sum.append(np.sum(CN_array * Eb_array))
            print(Eb_sum)
            Vr_max = []
            for Vr_dict in Vr:
                try:
                    Vr_max.append(max(Vr_dict.values()))
                except ValueError:
                    Vr_max.append(np.nan)
        except IndexError and ValueError:
            Eb_sum = np.nan
            Vr_max = np.nan
            print('no crystal features')
            pass
        mp_Eg = row["MP_Eg"]
        Eg = row["Eg"]
        mp_Ehull = row["MP_ehull"]
        Ev = row["vacancy_formation_energy"]
        bulk_type = row["bulk"]
        vacancy = row["group"]
        period = row["PT_oxide"]

        try:
            df_cf = pd.concat(
                [
                    df_cf,
                    pd.DataFrame(
                        {
                            "formula": formula,
                            "Ev": Ev,
                            "mp_Eg": mp_Eg,
                            "mp_ehull": mp_Ehull,
                            "Eg": Eg,
                            "bulk_type": bulk_type,
                            "vacancy": vacancy,
                            "space_group": space_group,
                            "Eb": Eb_sum,
                            "Vr_max": Vr_max,
                            "PT_oxide": period
                        }
                    )
                ]
            )
        except ValueError:
            print("df fail")
            pass
        df_cf.to_csv("wexler_nw_match.csv")

if __name__ == "__main__":
    main()
