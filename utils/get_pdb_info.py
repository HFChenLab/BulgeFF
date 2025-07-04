import logging

logging.basicConfig(filename='BulgeFix.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

"""
Extracts information from a PDB file and returns it based on the specified option.
    
Parameters:
pdb_name (str): The name of the PDB file to read.
option (str): The option specifying the type of information to return. 
            - 'full' returns a list of dictionaries with detailed atom information.
            - 'coord' returns a dictionary with atom coordinates.
    
Returns:
list or dict: 
    - If option is 'full': A list of dictionaries where each dictionary contains detailed information for each atom.
    - If option is 'coord': A dictionary with atom coordinates keyed by a tuple (residue_name, residue_number, atom_name).
"""

def get_pdb_info(pdb_name, option):
    pdb_info = []
    coords = {}
    
    try:
        with open(pdb_name, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith('ATOM'):
                    atom_number = int(line[6:11])
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    chain_id = line[21]
                    residue_number = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    occupancy = float(line[54:60])
                    bfactor = float(line[60:66])
                    element = line[76:78].strip()

                    pdb_info.append({
                        'atom_number': atom_number,
                        'atom_name': atom_name,
                        'residue_name': residue_name,
                        'chain_id': chain_id,
                        'residue_number': residue_number,
                        'x': x,
                        'y': y,
                        'z': z,
                        'occupancy': occupancy,
                        'bfactor': bfactor,
                        'element': element
                    })
                    
                    coords[(residue_name, residue_number, atom_name)] = (x, y, z)

        logging.info(f"Read PDB file: {pdb_name}.")

    except FileNotFoundError:
        error_msg = f"PDB file '{pdb_name}' not found. Please cheak the file path."
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)
    except Exception as e:
        error_msg = f"An error occurred: {str(e)}"
        logging.error(error_msg)
        raise

    if option == "full":
        return pdb_info
    elif option == "coord":
        return coords
    else:
        error_msg = f"Unknown option: {option}. Valid options are 'full' and 'coord'."
        logging.error(error_msg)
        raise ValueError(error_msg)
