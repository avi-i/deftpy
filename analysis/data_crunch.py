import os
import shutil
import subprocess
import bz2
import sys
from pymatgen.core import Structure
import glob

def calculate_bandgap(outcar_file):

    # Extract start and end lines from the command line arguments

    # Command to find the line number of the last occurrence of "occupation" in the OUTCAR file
    start_line_command = f"bzgrep -n 'occupation' '{outcar_file}' | tail -1 | awk '{{print $1}}' | sed 's|:||g'"
    start_line_output = subprocess.check_output(start_line_command, shell=True)
    start_line = int(start_line_output.strip())

    # Command to find the value of NBANDS in the OUTCAR file
    bands_command = f"bzgrep 'NBANDS=' '{outcar_file}' | tail -1 | awk '{{print $15}}'"
    bands_output = subprocess.check_output(bands_command, shell=True)
    bands = int(bands_output.strip())

    # Calculate the end line
    end_line = start_line + bands

    # Open the bz2 file
    with bz2.open(outcar_file, 'rt') as file:
        # Initialize variables to store the last flag 1 value and the first flag 0 value
        last_flag_1_value = None
        first_flag_0_value = None

        # Skip lines until start_line
        for _ in range(start_line - 1):
            file.readline()

        # Iterate through lines within the specified range
        for line_number in range(start_line, end_line + 1):
            line = file.readline()
            columns = line.split()
            if len(columns) == 3:
                value = float(columns[1])
                flag = float(columns[2])

                # Check if flag changes from 1 to 0
                if flag == 0.0 and last_flag_1_value is not None:
                    first_flag_0_value = value
                    break  # Exit the loop after finding the first flag 0 value

                # Update last_flag_1_value with the latest flag 1 value
                if flag == 1.0:
                    last_flag_1_value = value

    # Return the calculated bandgap
    return first_flag_0_value - last_flag_1_value

def keep_last_x_values(lst, x):
    return lst[-x:]

def calculate_bandgap_kpts(outcar_file):
    # Extract start and end lines from the command line arguments

    # Command to find the line number of the last occurrence of "occupation" in the OUTCAR file
    start_lines_command = f"bzgrep -n 'occupation' '{outcar_file}' | awk '{{print $1}}' | sed 's|:||g'"
    output = subprocess.check_output(start_lines_command, shell=True)

    # Decode the output and split it into lines
    output_lines = output.decode("utf-8").split("\n")

    # Remove empty lines
    output_lines = [line.strip() for line in output_lines if line.strip()]

    # Convert the lines to integers and create a Python list
    start_lines = [int(line) for line in output_lines]

    # kpoints
    last_occu_command = f"bzgrep -n 'occupation' '{outcar_file}' | tail -1 | awk '{{print $1}}' | sed 's|:||g'"
    last_occu_output = subprocess.check_output(last_occu_command, shell=True)
    last_occu = int(last_occu_output.strip())
    kpt_line = last_occu - 1

    with bz2.open(outcar_file, 'rt') as file:
        file.seek(0)  # Reset file pointer to the beginning of the file
        for _ in range(kpt_line - 1):
            file.readline()

        # Read the line at the desired line number
        line = file.readline()

        # Find the first integer in the line
        for word in line.split():
            try:
                integer_value = int(word)
                kpts = integer_value
            except ValueError:
                pass

    # Command to find the value of NBANDS in the OUTCAR file
    bands_command = f"bzgrep 'NBANDS=' '{outcar_file}' | tail -1 | awk '{{print $15}}'"
    bands_output = subprocess.check_output(bands_command, shell=True)
    bands = int(bands_output.strip())

    # Calculate the end line
    gaps = []
    for start_line in start_lines:
        end_line = start_line + bands

        # Open the bz2 file
        with bz2.open(outcar_file, 'rt') as file:
            # Initialize variables to store the last flag 1 value and the first flag 0 value
            last_flag_1_value = None
            first_flag_0_value = None

            # Skip lines until start_line
            for _ in range(start_line - 1):
                file.readline()

            # Iterate through lines within the specified range
            for line_number in range(start_line, end_line + 1):
                line = file.readline()
                columns = line.split()
                if len(columns) == 3:
                    value = float(columns[1])
                    flag = float(columns[2])

                    # Check if flag changes from 1 to 0
                    if flag == 0.0 and last_flag_1_value is not None:
                        first_flag_0_value = value
                        break  # Exit the loop after finding the first flag 0 value

                    # Update last_flag_1_value with the latest flag 1 value
                    if flag == 1.0:
                        last_flag_1_value = value
        gap = first_flag_0_value - last_flag_1_value
        gaps.append(gap)
    # Return the calculated bandgap
    new_gaps = keep_last_x_values(gaps, (kpts * 2))
    return min(new_gaps)

with open('/data3/wexler/data_collection_py_rerun.csv', 'w') as f:
    f.write("PT_oxide,formula,space_group,E_hi_sym_bulk,n_bulk_hi_sym,n_bulk,E_rattle_bulk,rattle_txt,hi_sym_vac,"
            "hi_sym_vac_BD,E_hi_sym_vac,bandgap_hi_sym_vac,rattle_vac_site,rattle_vac_BD,E_rattle_vac,"
            "bandgap_rattle_vac,MP_ehull,MP_Eg,vac_rerun,E_vac_rerun,vac_rerun_Eg,\n")

    binary_directories = [directory for directory in glob.glob('binary_*') if os.path.isdir(directory)]
    print(binary_directories)
    for p in binary_directories:
        os.chdir(p)
        print("inside", p)
        PT = p.replace('binary_', '').replace('/', '')
        dir = p.replace('/', '')
        print(PT)
        all_entries = os.listdir('.')
        # Filter out only the directories from the list of entries
        metal_oxides = [entry for entry in all_entries if os.path.isdir(entry)]
        for s in metal_oxides:
            os.chdir(s)
            print("inside", s)
            formula = s.split('_')[0]
            space_group = s[s.find('_') + 1:].replace('/', '')

            hs = [d for d in os.listdir() if d.startswith('higher')]
            for higher_sym in hs:
                os.chdir(higher_sym)
                print(f"in {higher_sym}")
                rerun_dirs = sorted([d for d in os.listdir() if d.startswith('rerun')],
                                    key=lambda x: int(x.split('-')[-1]) if '-' in x else int(x.replace('rerun', '0')),
                                    reverse=True)
                if rerun_dirs:
                    rerun = rerun_dirs[0]
                    os.chdir(rerun)
                    print(f"in {rerun}")
                else:
                    print("No rerun directory found, checking current directory")
                    rerun = os.getcwd()  # Use the current directory
                    print(f"in {rerun}")
                    os.chdir(rerun)
                E_hi_sym_bulk = float(subprocess.check_output(['bzgrep', 'energy  without entropy=', 'OUTCAR.bz2']).decode().split()[-1])
                src_path = os.path.join(os.getcwd(), "POSCAR.bz2")
                dst_folder = "/data3/wexler/poscars_py_rerun"
                dst_filename = f"{PT}_{s}_{higher_sym}_POSCAR.bz2"
                dst_path = os.path.join(dst_folder, dst_filename)
                if os.path.exists(src_path) and os.path.exists(dst_folder):
                    shutil.copyfile(src_path, dst_path)
                else:
                    print("Source file or destination folder does not exist.")
                with bz2.open(filename=src_path, mode="rb") as p:
                    hi_sym_bulk_poscar = p.read().decode('utf-8')
                    hi_sym_bulk_structure = Structure.from_str(hi_sym_bulk_poscar, fmt="poscar")
                    n_atom_hi_sym = hi_sym_bulk_structure.num_sites
                vacancies = [d for d in os.listdir() if d.endswith('vacancies')]
                if vacancies:
                    os.chdir(vacancies[0])
                    vacancy = [d for d in os.listdir() if d.startswith('v_')]
                    for v in vacancy:
                        os.chdir(v)
                        hi_sym_vac = v.replace('v_O_', '').replace('/', '')
                        bond_distortions = [d for d in os.listdir() if os.path.isdir(d)]
                        for bond in bond_distortions:
                            os.chdir(bond)
                            hi_sym_vac_BD = bond.split('_')[2]
                            src_path = os.path.join(os.getcwd(), "POSCAR.bz2")
                            dst_folder = "/data3/wexler/poscars_py_rerun"
                            dst_filename = f"{PT}_{s}_{hi_sym_vac}_{hi_sym_vac_BD}_POSCAR.bz2"
                            dst_path = os.path.join(dst_folder, dst_filename)
                            if os.path.exists(src_path) and os.path.exists(dst_folder):
                                shutil.copyfile(src_path, dst_path)
                            else:
                                print("Source file or destination folder does not exist.")
                            with bz2.open(filename=src_path, mode="rb") as p2:
                                hi_sym_vac_poscar = p2.read().decode('utf-8')
                                hi_sym_vac_structure = Structure.from_str(hi_sym_vac_poscar, fmt="poscar")
                                n_atom_vac = hi_sym_vac_structure.num_sites
                                n_bulk_hi_sym = (n_atom_vac + 1) / n_atom_hi_sym
                            E_hi_sym_vac = float(subprocess.check_output(['bzgrep', 'energy  without entropy=', 'OUTCAR.bz2']).decode().split()[-1])
                            bandgap_hi_sym_vac = calculate_bandgap("OUTCAR.bz2")
                            print("hi_sym vacancies are in rerun")
                            hs_vac_reruns = sorted([d for d in os.listdir() if d.startswith('rerun')],
                                key=lambda x: int(x.split('-')[-1]) if '-' in x else int(x.replace('rerun', '0')),
                                reverse=True)
                            if hs_vac_reruns:
                                hs_vac_rerun = hs_vac_reruns[0]
                                os.chdir(hs_vac_rerun)
                                E_vac_rerun = float(subprocess.check_output(
                                    ['bzgrep', 'energy  without entropy=', 'OUTCAR.bz2']).decode().split()[-1])
                                vac_rerun_Eg = calculate_bandgap_kpts("OUTCAR.bz2")
                                os.chdir("../")
                                f.write(
                                    f"{PT},{formula},{space_group},{E_hi_sym_bulk},{n_bulk_hi_sym},,,,{hi_sym_vac},{hi_sym_vac_BD},{E_hi_sym_vac},{bandgap_hi_sym_vac},,,,,,,True,{E_vac_rerun},{vac_rerun_Eg},\n")
                            f.write(f"{PT},{formula},{space_group},{E_hi_sym_bulk},{n_bulk_hi_sym},,,,{hi_sym_vac},{hi_sym_vac_BD},{E_hi_sym_vac},{bandgap_hi_sym_vac},,,,,,,,,,\n")

                            os.chdir('../')
                        os.chdir('../')
                    os.chdir('../')
                else:
                    print("hi_sym doesn't have vacancies")
                    f.write(f"{PT},{formula},{space_group},{E_hi_sym_bulk},,,,,,,,,,,,,,,,\n")
                os.chdir(f'/data3/wexler/cfm/{dir}/{s}')

            os.chdir('rattle')
            print(f"in {dir} {s} rattle")
            rattle_txt = subprocess.check_output(['bzcat', 'rattle.txt.bz2']).decode()
            MP_Eg = subprocess.check_output(['bzcat', 'band_gap.txt.bz2']).decode().strip()
            MP_ehull = subprocess.check_output(['bzcat', 'energy_above_hull.txt.bz2']).decode().strip()
            rerun_dirs = sorted([d for d in os.listdir() if d.startswith('rerun')],
                                key=lambda x: int(x.split('-')[-1]) if '-' in x else int(x.replace('rerun', '0')),
                                reverse=True)
            if rerun_dirs:
                rerun = rerun_dirs[0]
                os.chdir(rerun)
                print(f"in {rerun}")
            else:
                print("No rerun directory found, checking current directory")
                rerun = os.getcwd()  # Use the current directory
                print(f"in {rerun}")
                os.chdir(rerun)
            E_rattle_bulk = float(subprocess.check_output(['bzgrep', 'energy  without entropy', 'OUTCAR.bz2']).decode().split()[-1])
            src_path = os.path.join(os.getcwd(), "POSCAR.bz2")
            dst_folder = "/data3/wexler/poscars_py_rerun"
            dst_filename = f"{PT}_{s}_POSCAR.bz2"
            dst_path = os.path.join(dst_folder, dst_filename)
            if os.path.exists(src_path) and os.path.exists(dst_folder):
                shutil.copyfile(src_path, dst_path)
            else:
                print("Source file or destination folder does not exist.")
            with bz2.open(filename=src_path, mode="rb") as p3:
                bulk_poscar = p3.read().decode('utf-8')
                bulk_structure = Structure.from_str(bulk_poscar, fmt="poscar")
                n_atom_bulk = bulk_structure.num_sites

            vacancies = [d for d in os.listdir() if d.endswith('vacancies')]
            if vacancies:
                os.chdir(vacancies[0])
                rattle_vacancies = [d for d in os.listdir() if d.startswith('v_')]
                for v in rattle_vacancies:
                    os.chdir(v)
                    rattle_vac_site = v.replace('v_O_', '').replace('/', '')

                    BD = [d for d in os.listdir() if os.path.isdir(d)]
                    for bond in BD:
                        os.chdir(bond)
                        rattle_vac_BD = bond.split('_')[2]
                        src_path = os.path.join(os.getcwd(), "POSCAR.bz2")
                        dst_folder = "/data3/wexler/poscars_py_rerun"
                        dst_filename = f"{PT}_{s}_{rattle_vac_site}_{rattle_vac_BD}_POSCAR.bz2"
                        dst_path = os.path.join(dst_folder, dst_filename)
                        if os.path.exists(src_path) and os.path.exists(dst_folder):
                            shutil.copyfile(src_path, dst_path)
                        else:
                            print("Source file or destination folder does not exist.")
                        E_rattle_vac = float(subprocess.check_output(['bzgrep', 'energy  without entropy=', 'OUTCAR.bz2']).decode().split()[-1])
                        bandgap_rattle_vac = calculate_bandgap("OUTCAR.bz2")
                        with bz2.open(filename=src_path, mode="rb") as p4:
                            rattle_vac_poscar = p4.read().decode('utf-8')
                            rattle_vac_structure = Structure.from_str(rattle_vac_poscar, fmt="poscar")
                            n_atom_rattle_vac = rattle_vac_structure.num_sites
                            n_bulk = (n_atom_rattle_vac + 1) / n_atom_bulk
                        print("rattle rerun has vacancies")
                        vac_reruns = sorted([d for d in os.listdir() if d.startswith('rerun')],
                                key=lambda x: int(x.split('-')[-1]) if '-' in x else int(x.replace('rerun', '0')),
                                reverse=True)
                        if vac_reruns:
                            vac_rerun = vac_reruns[0]
                            os.chdir(vac_rerun)
                            E_vac_rerun = float(subprocess.check_output(
                                ['bzgrep', 'energy  without entropy=', 'OUTCAR.bz2']).decode().split()[-1])
                            rat_vac_rerun_Eg = calculate_bandgap_kpts("OUTCAR.bz2")
                            os.chdir("../")
                            f.write(
                                f"{PT},{formula},{space_group},,,{n_bulk},{E_rattle_bulk},,,,,,{rattle_vac_site},{rattle_vac_BD},{E_rattle_vac},{bandgap_rattle_vac},{MP_ehull},{MP_Eg},True,{E_vac_rerun},{rat_vac_rerun_Eg},\n")
                        f.write(f"{PT},{formula},{space_group},,,{n_bulk},{E_rattle_bulk},,,,,,{rattle_vac_site},{rattle_vac_BD},{E_rattle_vac},{bandgap_rattle_vac},{MP_ehull},{MP_Eg},,,,\n")

                        os.chdir('../')
                    os.chdir('../')
                os.chdir('../')
            else:
                print("rattle rerun does not have vacancies")
                f.write(f"{PT},{formula},{space_group},,,,{E_rattle_bulk},,,,,,,,,,{MP_ehull},{MP_Eg},,,,\n")

            os.chdir(f'/data3/wexler/cfm/{dir}')

        os.chdir(f'/data3/wexler/cfm/')