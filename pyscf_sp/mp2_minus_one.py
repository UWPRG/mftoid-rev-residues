import os
from pyscf import gto, mp

# Path to the folder containing the xyz files
folder_path = 'minus_xyz/'

# Path to the output text file
output_file_path = 'energies.txt'


def get_structure(xyz_file):
    """
    Reads an XYZ file and returns the atomic coordinates.
    
    Parameters:
    xyz_file (str): The path to the XYZ file.
    
    Returns:
    list: A list of atomic coordinates.
    """
    lines = []
    with open(xyz_file) as xyz:
        for line in xyz:
            lines.append(line)

    return lines[2:]  # Skip the first two lines (header and atom count)


files = os.listdir(folder_path)

with open(output_file_path, 'w') as output_file:

    for file_name in files:
        
        xyz = os.path.join(folder_path, file_name)

        # Create a PySCF molecule object from the XYZ file
        mol = gto.M(
                atom=get_structure(xyz),  # Get atomic structure from the XYZ file
                basis="6-31G**",          # Basis set
                verbose=4,                # Verbosity level
                output="PySCF.log",       # Log file
                spin=0,                   # Spin multiplicity
                unit="Angstrom",          # Unit of distance
                symmetry=False,           # Symmetry (not used)
                charge=-1)                 # Molecular charge

        # Perform Hartree-Fock calculation
        hf = mol.RHF().run()

        # Perform MP2 calculation based on the Hartree-Fock result
        mp2 = mp.MP2(hf)
        mp2.kernel()

        # Extract the total MP2 energy
        energy = mp2.e_tot
        
        output_file.write(f'{file_name}\t{energy}\n')
        
        output_file.flush()
        os.remove('PySCF.log')
        
        print("Done")

