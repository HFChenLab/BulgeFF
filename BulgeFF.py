import sys
import os
import argparse
import warnings
from Bio import BiopythonDeprecationWarning

warnings.simplefilter('ignore', BiopythonDeprecationWarning)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils import (
    get_atom_id,
    get_trinucleotides,
    get_prototype,
    get_function,
    sugar_type,
    write_plumed
)

def main(bulge_pdb_name, md_pdb_name, bulge_resi_list, bulge_res_list, prototype_path, func_file, plumed_file):
    eta_list = []
    theta_list = []
    eta_func_list = []
    theta_func_list = []

    for bulge_resi, bulge_res in zip(bulge_resi_list, bulge_res_list):
        reference = get_trinucleotides(bulge_pdb_name, bulge_resi, bulge_res)
        eta, theta = get_atom_id(md_pdb_name, bulge_resi, bulge_res)
        sugar = sugar_type(bulge_pdb_name, bulge_resi, bulge_res)
        model, model_type = get_prototype(sugar, reference, prototype_path)

        if model_type is None:
            print(model)
            continue

        eta_func, theta_func = get_function(model_type, model, func_file)

        eta_list.append(eta)
        theta_list.append(theta)
        eta_func_list.append(eta_func)
        theta_func_list.append(theta_func)

    write_plumed(plumed_file, eta_list, theta_list, eta_func_list, theta_func_list, bulge_res_list, bulge_resi_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process bulge RNA data and write PLUMED input.",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True
    )

    parser.add_argument("--bulge_pdb", type=str, default="2jym.pdb", help="Name of PDB file. Default: 2jym.pdb")
    parser.add_argument("--md_pdb", type=str, default="reference.pdb", help="Name of PDB file for MD. Default: reference.pdb")
    parser.add_argument("--bulge_name", type=str, nargs='+', default=["G"], help="List of bulge residue names. Default: G")
    parser.add_argument("--bulge_id", type=int, nargs='+', default=[6], help="List of bulge residue IDs. Default: 6")

    args = parser.parse_args()

    if len(args.bulge_name) != len(args.bulge_id):
        print("Error: The number of bulge residues and residue IDs must match.")
        sys.exit(1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, 'function')

    prototype_path = os.path.join(script_dir, 'prototype_db')
    func_file = os.path.join(data_dir, 'fix_function.txt')
    plumed_file = os.path.join(os.getcwd(), 'plumed.dat')

    main(args.bulge_pdb, args.md_pdb, args.bulge_id, args.bulge_name, prototype_path, func_file, plumed_file)