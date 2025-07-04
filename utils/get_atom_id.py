from utils.get_pdb_info import get_pdb_info
import logging

logging.basicConfig(filename='BulgeFix.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

"""
Retrieves the atom IDs for a specified bulge residue and its neighboring residues from a PDB file.

Parameters:
pdb_name (str): The name of the PDB file to be read.
bulge_resi (int): The residue number of the bulge residue.
bulge_res (str): The name of the bulge residue.

Returns:
tuple: A tuple containing two strings:
    - `eta`: A string with atom IDs for the eta dihedral angle in the format "atom_id1,atom_id2,atom_id3,atom_id4".
    - `theta`: A string with atom IDs for the theta dihedral angle in the format "atom_id1,atom_id2,atom_id3,atom_id4".
"""

def get_atom_id(pdb_name, bulge_resi, bulge_res):

    pdb_info = get_pdb_info(pdb_name, option="full")

    logging.info(f"Processing: bulge residue: {bulge_res}, bulge residue id: {bulge_resi}")

    atom_id = {}
    chain_id = None
    found_residue = False

    for atom in pdb_info:
        if atom['residue_number'] == bulge_resi and atom['residue_name'] == bulge_res:
            chain_id = atom['chain_id']
            found_residue = True
            break
    
    if not found_residue:
        error_msg = f"Residue {bulge_res} with ID {bulge_resi} not found in any chain."
        logging.error(error_msg)
        raise ValueError(error_msg)

    for atom in pdb_info:
        if atom['chain_id'] == chain_id:
            if atom['residue_number'] == bulge_resi and atom['residue_name'] == bulge_res:
                if atom['atom_name'] == "P":
                    atom_id['P_i'] = (atom['atom_number'], atom['residue_number'])
                elif atom['atom_name'] == "C4'":
                    atom_id['C4_i'] = (atom['atom_number'], atom['residue_number'])
            elif atom['residue_number'] == bulge_resi - 1:
                if atom['atom_name'] == "C4'":
                    atom_id['C4_i-1'] = (atom['atom_number'], atom['residue_number'])
            elif atom['residue_number'] == bulge_resi + 1:
                if atom['atom_name'] == "P":
                    atom_id['P_i+1'] = (atom['atom_number'], atom['residue_number'])
                elif atom['atom_name'] == "C4'":
                    atom_id['C4_i+1'] = (atom['atom_number'], atom['residue_number'])

    required_atoms = {
        'C4_i-1': f"C4'_{bulge_resi-1}",
        'P_i': f"P_{bulge_resi}",
        'C4_i': f"C4'_{bulge_resi}",
        'P_i+1': f"P_{bulge_resi+1}",
        'C4_i+1': f"C4'_{bulge_resi+1}"
    }
    
    missing_atoms = [required_atoms[atom] for atom in required_atoms if atom_id.get(atom) is None]
    if missing_atoms:
        error_msg = f"Missing atom(s) for eta/theta: {', '.join(missing_atoms)}"
        logging.error(error_msg)
        raise ValueError(error_msg)
    
    eta = f"{atom_id['C4_i-1'][0]},{atom_id['P_i'][0]},{atom_id['C4_i'][0]},{atom_id['P_i+1'][0]}"
    theta = f"{atom_id['P_i'][0]},{atom_id['C4_i'][0]},{atom_id['P_i+1'][0]},{atom_id['C4_i+1'][0]}"
    
    logging.info(f"Atom id of eta: {eta}")
    logging.info(f"Atom id of theta: {theta}")
    
    return eta, theta
