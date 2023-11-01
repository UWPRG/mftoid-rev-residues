import mdtraj as md
import numpy as np
import networkx as nx
import Bio.PDB
import os
import argparse
import subprocess
import shutil

# from utils import *

filenames = {
    'A': 'ALAp',
    'R': 'ARGp',
    'N': 'ASNp',
    'D': 'ASPd',
    'C': 'CYSd',
    'E': 'GLUp',
    'Q': 'GLNp',
    'G': 'GLYp',
    'H': 'HSDp',
    'J': 'HSEp',
    'B': 'HSPp',
    'I': 'ILEp',
    'L': 'LEUp',
    'K': 'LYSp',
    'M': 'METd',
    'F': 'PHEp',
    'P': 'PROd',
    'S': 'SERp',
    'T': 'THRp',
    'W': 'TRPp',
    'Y': 'TYRp',
    'V': 'VALp'
}


parser = argparse.ArgumentParser(
    description="Molecule Inserter for peptoids"
)
parser.add_argument(
    "--seq",
    type=str,
    default=None,
    metavar="seq",
    help="box size",
)
parser.add_argument(
    "--bs",
    type=str,
    default=None,
    metavar="bs",
    help="box size",
)
args = parser.parse_args()

box_size = args.bs
if box_size is None:
    box_size = 3
else:
    try:
        box_size = float(box_size)
    except:
        raise ValueError("Box size must be a number")
if args.seq is None:
    raise ValueError("Sequence Entry Needed")
    
peptoid_size = len(args.seq)


nterm = md.load_pdb("ACEp_f.pdb", standard_names=False)
cterm = md.load_pdb("NMEp_f.pdb", standard_names=False)
os.mkdir(os.path.join(os.getcwd(), "med"))
command = "obabel ACEp_f.pdb -O med/ace_h.pdb -h"
subprocess.run(command.split(), stdout=subprocess.PIPE)
nh = md.load("med/ace_h.pdb")
prev_coord = nh.xyz[0, -1]

full_structure = md.Trajectory(nterm.xyz, nterm.topology)

for letter in args.seq:
    if letter not in filenames.keys():
        err_string = str(letter) + " is not a valid amino acid code."
        raise ValueError(err_string)
    # Read the structure PDB file
    filename = "residue_pdb/" + filenames[letter] + ".pdb"
    structure = md.load(filename, standard_names=False)
    #Read in the proper orientation of alanine
    orientation = md.load("orientation.pdb")
    or_indices = [0, 6, 8, 9] #indices of the nitrogen, alpha carbon, carbonyl carbon, and carbonyl oxygen in ALAp template file.
    mid_indices = structure.topology.select("name NL or name CA or name CLP or name OL")
    #Superimpose your amino acid onto alanine
    structure.superpose(orientation, atom_indices=mid_indices, ref_atom_indices=or_indices)

    #Save as new pdb
    structure.save_pdb('med/rotated.pdb')
    #Add hydrogens
    command = "obabel med/rotated.pdb -O med/rotated_h.pdb -h"
    subprocess.run(command.split(), stdout=subprocess.PIPE)

    #Collect these hydrogens' positions
    rot = md.load("med/rotated_h.pdb")
    n_atoms = rot.n_atoms
    nterm_coord = rot.xyz[0, n_atoms - 2]
    translate = prev_coord - nterm_coord
    structure.xyz[0] += translate
    rot.xyz[0] += translate
    prev_coord = rot.xyz[0, n_atoms - 1]

    full_structure = full_structure.stack(structure)





#Find where the n-terminus and the c-terminus are bonding
# n_index = nterm.topology.select("name CLP")[0]
# n_coord = nterm.xyz[0, n_index]
c_index = cterm.topology.select("name NL")[0]
c_coord = cterm.xyz[0, c_index]

#Find where to place n-terminus and c-terminus according to where the bonding atoms currently are and where they should be
c_translate = prev_coord - c_coord

cterm.xyz += c_translate
full_structure = full_structure.stack(cterm)

#load the pdb without hydrogens again
#stack the amino acid onto the n-terminus
top = full_structure.topology
#make new bonds




# Iterate over all residues in the topology

        
# cterm_nl_index = new_traj.topology.select("name NL")[1]


table, bonds = top.to_dataframe()
for i in range(1, len(table['chainID'])):
    mask = table['chainID'] == i
    table.loc[mask, 'resSeq'] = i
table[:]['chainID'] = 0
table[:]['serial'] = 1 + np.arange(full_structure.n_atoms)


# table[:]['resName'] += "p"

new_top = md.Topology.from_dataframe(table, bonds)

aa_nl_indices = new_top.select("name NL")
aa_clp_indices = new_top.select("name CLP")

atoms = [atom for atom in new_top.atoms]

for i in range(len(aa_nl_indices)):
    new_top.add_bond(atoms[aa_nl_indices[i]], atoms[aa_clp_indices[i]])
full_structure.topology = new_top

full_structure.save_pdb("combined.pdb")


# command = "obabel combined.pdb -O output.pdb"
# subprocess.run(command.split(), stdout=subprocess.PIPE)
with open("combined.pdb") as f:
    lines = f.readlines()

with open("output.pdb", "w") as f:
    for line in lines:
        if line.startswith("ATOM"):
            if 'PRO' in line or 'GLU' in line or 'ASP' in line or 'CYS' in line or 'MET' in line:
                l = line[:20] + "d    " + line[25:]
            else:
                l = line[:20] + "p    " + line[25:]
        else:
            l = line
        f.write(l)


shutil.rmtree(os.path.join(os.getcwd(), "med"))
