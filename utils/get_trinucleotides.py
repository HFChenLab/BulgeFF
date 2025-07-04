from utils.get_pdb_info import get_pdb_info
from utils.atom_to_pdb import atom_to_pdb
import logging

logging.basicConfig(filename='BulgeFix.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

"""
Extracts trinucleotide information from a PDB file based on the specified bulge residue and its neighboring residues.

Parameters:
- pdb_name (str): The name of the PDB file to be read.
- bulge_resi (int): The residue number of the bulge.
- bulge_res (str): The residue name of the bulge.

Returns:
- trinucleotides (str): A PDB-format string representing the atoms of the trinucleotide structure including the bulge and its neighboring residues.
"""

def get_trinucleotides(pdb_name, bulge_resi, bulge_res):
    try:
        pdb_info = get_pdb_info(pdb_name, option="full")
        logging.info(f"Processing trinucleotides with bulge residue: {bulge_resi}-{bulge_res}")

        target_residues = {bulge_resi - 1, bulge_resi, bulge_resi + 1}
        atom_info = []
        residue_names = {}
        found_residues = set()

        for atom in pdb_info:
            if (atom['residue_number'] in target_residues and
                (atom['residue_number'] != bulge_resi or atom['residue_name'] == bulge_res)):
                atom_info.append(atom)
                found_residues.add(atom['residue_number'])
                residue_names[atom['residue_number']] = atom['residue_name']

        missing_residues = target_residues - found_residues
        if missing_residues:
            error_msg = f"Missing residues in target set: {missing_residues}"
            logging.error(error_msg)
            raise ValueError(error_msg)

        found_residues_info = {res_num: residue_names[res_num] for res_num in found_residues}
        logging.info(f"Found trinucleotide atoms for residues: {found_residues_info}")

        trinucleotides = atom_to_pdb(atom_info)

    except Exception as e:
        logging.error(f"An error occurred while processing trinucleotides: {str(e)}")
        raise

    return trinucleotides
