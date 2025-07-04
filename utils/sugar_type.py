from utils.get_pdb_info import get_pdb_info
import math
import numpy as np
import logging


logging.basicConfig(filename='BulgeFix.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

"""
Determines the sugar type (C2'-endo, C3'-endo, or Others) based on the phase angle calculated from PDB coordinates.

Parameters:
- pdb_name (str): The name of the PDB file containing the structure.
- bulge_resi (int): The residue identifier of the bulge.
- bulge_res (str): The residue name of the bulge.

Returns:
- sugar_type (str): The determined sugar type based on the phase angle calculation.
"""

def calculate_dihedral_angle(p1, p2, p3, p4):
    b1 = p1 - p2
    b2 = p3 - p2
    b3 = p4 - p3

    b2 /= np.linalg.norm(b2)

    v = b1 - np.dot(b1, b2) * b2
    w = b3 - np.dot(b3, b2) * b2

    x = np.dot(v, w)
    y = np.dot(np.cross(b2, v), w)

    return math.atan2(y, x)

def sugar_type(pdb_name, bulge_resi, bulge_res):
    try:
        logging.info(f"Calculating phase angle for residue {bulge_res}-{bulge_resi} in file {pdb_name}")

        coords = get_pdb_info(pdb_name, option="coord")

        C4 = np.array(coords.get((bulge_res, bulge_resi, "C4'")))
        O4 = np.array(coords.get((bulge_res, bulge_resi, "O4'")))
        C1 = np.array(coords.get((bulge_res, bulge_resi, "C1'")))
        C2 = np.array(coords.get((bulge_res, bulge_resi, "C2'")))
        C3 = np.array(coords.get((bulge_res, bulge_resi, "C3'")))
        
        missing_atoms = []
        for atom_name, coords in [("C4'", C4), ("O4'", O4), ("C1'", C1), ("C2'", C2), ("C3'", C3)]:
            if None in coords:
                missing_atoms.append(atom_name)
        
        if missing_atoms:
            error_msg = (f"Missing atoms for residue {bulge_res}-{bulge_resi} when calculating phase angle: "
                         f"{', '.join(missing_atoms)}")
            raise ValueError(error_msg)

        v0 = calculate_dihedral_angle(C4, O4, C1, C2)
        v1 = calculate_dihedral_angle(O4, C1, C2, C3)
        v2 = calculate_dihedral_angle(C1, C2, C3, C4)
        v3 = calculate_dihedral_angle(C2, C3, C4, O4)
        v4 = calculate_dihedral_angle(C3, C4, O4, C1)

        v = [v2, v3, v4, v0, v1]

        A = 0
        B = 0
        for j in range(1, 6):
            t = 0.8 * math.pi * (j - 1)
            A += v[j-1] * math.cos(t)
            B += v[j-1] * math.sin(t)

        A *= 0.4
        B *= -0.4

        tm = math.sqrt(A * A + B * B)

        c = A / tm
        s = B / tm

        phase = math.atan2(s, c) * 180 / math.pi
        if phase < 0:
            phase += 360

        if (phase >= 0 and phase <= 54) or (phase >= 342 and phase <= 360):
            sugar_type = "C3'-endo"
        elif phase >= 126 and phase <= 198:
            sugar_type = "C2'-endo"
        else:
            sugar_type = "Others"

        logging.info(f"Calculated phase angle: {phase} degrees, sugar type: {sugar_type}")
        return sugar_type

    except ValueError as e:
        logging.error(f"An error occurred while calculating the phase angle: {str(e)}")
        raise
    except Exception as e:
        logging.error(f"An unexpected error occurred: {str(e)}")
        raise