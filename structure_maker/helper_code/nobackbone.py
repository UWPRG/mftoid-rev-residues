import mdtraj as md
import numpy as np
import sys

def remove_backbone(filename):
    mol = md.load(filename)
    correct_indices = np.array(list(mol.topology.select('not name NL and not name CA and not name CLP and not name OL and not name HA1 and not name HA2')))
    table, bonds = mol.topology.to_dataframe()
    onswitch = np.zeros(len(table), dtype=np.bool)
    onswitch[correct_indices] = 1
    new_table = table[onswitch].reset_index()
    new_table['index'] = np.arange(len(new_table))
    new_table['serial'] = 1 + np.arange(len(new_table))
    new_top = md.Topology.from_dataframe(new_table)
    new_traj = md.Trajectory(mol.xyz[0, onswitch], new_top)
    new_traj.save_pdb(filename)
    
remove_backbone(str(sys.argv[1]))