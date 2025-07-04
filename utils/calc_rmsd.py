import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from io import StringIO
import logging
import warnings

logging.basicConfig(filename='BulgeFix.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

"""
Calculates the Root Mean Square Deviation (RMSD) between a reference structure and a prototype structure.
    
Parameters:
reference (str): The reference PDB structure as a string.
prototype (str): The prototype PDB structure as a string.
    
Returns:
float: The RMSD value between the reference and prototype structures.
"""

def calc_rmsd(reference, prototype):
    try:
        logging.info("Starting RMSD calculation.")
        reference_pdb = StringIO(reference)
        prototype_pdb = StringIO(prototype)

        u1 = mda.Universe(reference_pdb, format='pdb')
        u2 = mda.Universe(prototype_pdb, format='pdb')

        ref_atoms = u1.select_atoms("name C4' P")
        prototype_atoms = u2.select_atoms("name C4' P")

        logging.info(f"Selected atoms for alignment: {len(ref_atoms)} atoms in reference, {len(prototype_atoms)} atoms in prototype.")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            aligner = align.AlignTraj(u1, u2, select="name C4' P", in_memory=True)
            aligner.run()

            rmsd_value = rmsd(prototype_atoms.positions, ref_atoms.positions, superposition=True)

        logging.info(f"RMSD calculation completed: {rmsd_value:.3f} Å")

        return rmsd_value

    except Exception as e:
        logging.error(f"An error occurred during RMSD calculation: {str(e)}")
        raise
