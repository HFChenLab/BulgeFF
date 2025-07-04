import logging

logging.basicConfig(filename='BulgeFix.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

"""
Reads a function file to retrieve the eta and theta functions based on the provided sugar type and prototype.

Parameters:
- sugar_type (str): The type of sugar for which the function is to be retrieved.
- prototype (str): The prototype name associated with the function.
- func_file (str): The path to the file containing the functions.

Returns:
- tuple: A tuple containing two elements:
    - `eta_function` (str or None): The eta function corresponding to the provided sugar type and prototype, or None if not found.
    - `theta_function` (str or None): The theta function corresponding to the provided sugar type and prototype, or None if not found.
"""

def get_function(sugar_type, prototype, func_file):
    eta_function = None
    theta_function = None

    try:
        with open(func_file, 'r') as file:
            lines = file.readlines()

            for line in lines[1:]:
                parts = line.strip().split()
                if len(parts) == 4:
                    st, pt, dihedral, function = parts
                    if st == sugar_type and pt == prototype:
                        if dihedral == 'eta':
                            eta_function = function
                            logging.info(f"Found eta function for {sugar_type}, {prototype}: {function}")
                        elif dihedral == 'theta':
                            theta_function = function
                            logging.info(f"Found theta function for {sugar_type}, {prototype}: {function}")

            if eta_function is None:
                logging.warning(f"No eta function found for {sugar_type}, {prototype}")
            if theta_function is None:
                logging.warning(f"No theta function found for {sugar_type}, {prototype}")

    except FileNotFoundError:
        logging.error(f"Function file not found: {func_file}")
        raise
    except Exception as e:
        logging.error(f"Error reading function file {func_file}: {e}")
        raise

    return eta_function, theta_function