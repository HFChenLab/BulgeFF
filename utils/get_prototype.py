from Bio import PDB
from io import StringIO
from utils.calc_rmsd import calc_rmsd
import logging

logging.basicConfig(filename='BulgeFix.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

"""
Finds and selects the best prototype model based on RMSD value from the provided prototypes.

Parameters:
- sugar (str): The type of sugar (e.g., "Others", "C2'-endo", "C3'-endo") to determine which prototype files to load.
- reference (str): The PDB-format string of the reference structure for RMSD calculation.
- prototype_path (str): The directory path where prototype PDB files are located.

Returns:
- model (str): The name of the selected model based on the lowest RMSD value.
- model_type (str): The type of the selected model (e.g., "C2'-endo", "C3'-endo", or the value of `sugar`).
"""

def get_prototype_models(pdb_name):
    try:
        with open(pdb_name, 'r') as file:
            pdb_lines = file.readlines()
    except Exception as e:
        logging.error(f"Error reading PDB file {pdb_name}: {e}")
        raise

    model_indices = []
    models = []
    current_model = []
    for line in pdb_lines:
        if line.startswith("MODEL"):
            if current_model:
                models.append(''.join(current_model))
                current_model = []
            model_indices.append(int(line.split()[1]))
        current_model.append(line)
    
    if current_model:
        models.append(''.join(current_model))
    
    logging.info(f"Successfully parsed {len(models)} models from {pdb_name}")
    return models, model_indices


def get_prototype(sugar, reference, prototype_path):
    threshold = {"C2'-endo": 1.3, "C3'-endo": 1.2}
    
    try:
        if sugar == "Others":
            logging.info("Sugar type is 'Others', loading both C2'-endo and C3'-endo prototypes.")
            prototypes_C2, indices_C2 = get_prototype_models(f"{prototype_path}/C2'-endo.pdb")
            prototypes_C3, indices_C3 = get_prototype_models(f"{prototype_path}/C3'-endo.pdb")
            prototypes = prototypes_C2 + prototypes_C3
            indices = indices_C2 + indices_C3
        else:
            prototype_file = f"{prototype_path}/{sugar}.pdb"
            prototypes, indices = get_prototype_models(prototype_file)
    except Exception as e:
        logging.error(f"Error loading prototypes for sugar type {sugar}: {e}")
        raise
    
    rmsd_values = [calc_rmsd(reference, prototype) for prototype in prototypes]
    lowest_model_index = min(enumerate(rmsd_values), key=lambda x: x[1])[0]
    lowest_rmsd_value = rmsd_values[lowest_model_index]
    
    actual_model_index = indices[lowest_model_index]
    
    if sugar == "Others":
        model_type = "C2'-endo" if lowest_model_index < len(prototypes_C2) else "C3'-endo"
    else:
        model_type = sugar

    logging.info(f"Lowest RMSD value for {model_type}: {lowest_rmsd_value:.3f} Å")

    if model_type in threshold and lowest_rmsd_value > threshold[model_type]:
        warning_message = f"Error: No suitable prototype found as the RMSD value exceeds the threshold for {model_type}. RMSD: {lowest_rmsd_value:.3f} Å, Threshold: {threshold[model_type]} Å"
        logging.warning(warning_message)
        return warning_message, None
    
    model = f"model{actual_model_index}"
    logging.info(f"Selected model: {model} for sugar type: {model_type}")
    return model, model_type
