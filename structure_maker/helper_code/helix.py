import numpy as np
import mdtraj as md
import Bio.PDB
import os
import argparse
import subprocess
import shutil

from Bio.PDB.vectors import Vector, rotaxis, calc_dihedral, rotmat

from rdkit import Chem
from rdkit.Chem import AllChem

# from utils import *
kT = 50
filenames = {
    'A': 'ALAp',
    'R': 'ARGp',
    'N': 'ASNp',
    'D': 'ASPd',
    'C': 'CYSd',
    'E': 'GLUd',
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
num_iters = {
    'A': 0,
    'R': 100,
    'N': 100,
    'D': 10,
    'C': 10,
    'E': 10,
    'Q': 100,
    'G': 0,
    'H': 100,
    'J': 100,
    'B': 100,
    'I': 100,
    'L': 100,
    'K': 100,
    'M': 100,
    'F': 500,
    'P': 0,
    'S': 10,
    'T': 100,
    'W': 500,
    'Y': 500,
    'V': 100
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

args = parser.parse_args()

if args.seq is None:
    raise ValueError("Sequence Entry Needed")

backbone = md.load("backbone_pdb/backbone.pdb")
template = md.load("backbone_pdb/backbone_hn.pdb")


full_structure = md.Trajectory(np.zeros((1, 0, 3)), md.Topology())
nres = backbone.topology.n_residues
first_indices = backbone.topology.select("resid 0 or resid 1")
last_indices = backbone.topology.select("resid 0 or resid 1")
first_mean = np.mean(backbone.xyz[0, first_indices], axis=0)
last_mean = np.mean(backbone.xyz[0, last_indices], axis=0)
linear_chain = last_mean - first_mean
backbone_indices = []
sidechain_indices = []
cur_index = 0
i = 1
for letter in args.seq:
    backbone_indices.append(np.arange(cur_index, cur_index + 4))
    if letter not in filenames.keys():
        err_string = str(letter) + " is not a valid amino acid code."
        raise ValueError(err_string)
    # Read the structure PDB file
    filename = "noh_residue_pdb/" + filenames[letter] + ".pdb"
    cur_index += 4
    if letter == "P":
        # I hate proline!
        orientation = md.load('PROd.pdb', standard_names=False)
        template_backbone_indices = template.topology.select('resid ' + str(i-1) + " and backbone")
        pro_backbone_indices = [3, 1, 0, 2]
        pro_sidechain_indices = [4, 5, 8]
        orientation.superpose(template, atom_indices=pro_backbone_indices, ref_atom_indices=template_backbone_indices)
        structure = md.load(filename, standard_names=False)
        cur_index += structure.topology.n_atoms
        structure.xyz[0] = orientation.xyz[0, pro_sidechain_indices]
        template_ca_index = template_backbone_indices = template.topology.select('resid ' + str(i-1) + " and name CA")
        cb_index = structure.topology.select('name CB')
        ca_cb_bond = template.xyz[0, template_ca_index] - structure.xyz[0, cb_index]
        structure.xyz[0, cb_index] += 0.25 * ca_cb_bond
        t, b = structure.topology.to_dataframe()
        t['resSeq'] = i
        structure = md.Trajectory(structure.xyz, md.Topology.from_dataframe(t, b))
        backbone = backbone.stack(structure)
        table, bonds = backbone.topology.to_dataframe()
        resname = list(table.loc[table['chainID'] == 1, 'resName'])[0]
        table['chainID'] = 0
        table.loc[table['resSeq'] > nres, 'resSeq'] = i
        table.loc[table['resSeq'] == i, 'resName'] = resname
        table.loc[:, 'serial'] = np.arange(len(table)) + 1
        top = md.Topology.from_dataframe(table, bonds)
        t2, b2 = top.to_dataframe()
        proper_indices = t2['serial'].to_numpy() - 1
        backbone = md.Trajectory(backbone.xyz[0, proper_indices], top)
        top = backbone.topology
        bond_n = top.select("resid " + str(i-1) + " and name N")[0]
        bond_ca = top.select("resid " + str(i-1) + " and name CA")[0]
        bond_r = top.select("resSeq " + str(i) + " and resname " + resname + " and not name C and not name N and not name CA and not name O")[0]
        bond_c = top.select("resSeq " + str(i) + " and resname " + resname + " and not name C and not name N and not name CA and not name O")[0]
        atoms = [atom for atom in top.atoms]
        top.add_bond(atoms[np.where(proper_indices == bond_n)[0][0]], atoms[np.where(proper_indices == bond_r)[0][0]])
        top.add_bond(atoms[np.where(proper_indices == bond_ca)[0][0]], atoms[np.where(proper_indices == bond_c)[0][0]])
        backbone.topology = top
    elif letter == "G":
        sidechain_indices.append(np.array([], dtype=np.int32))
        table, bonds = backbone.topology.to_dataframe()
        table.loc[table['resSeq'] == i, 'resName'] = "GLY"
        backbone.topology = md.Topology.from_dataframe(table, bonds)
    else:
        structure = md.load(filename, standard_names=False)
        sidechain_indices.append(np.arange(structure.topology.n_atoms) + cur_index)
        cur_index += structure.topology.n_atoms
        cur_h = template.topology.select("resid " + str(i-1) + " and name H")[0]
        cur_n = template.topology.select("resid " + str(i-1) + " and name N")[0]
        if letter != "A":
            first_bond = template.xyz[0, cur_h] - template.xyz[0, cur_n]
            cur_orientation = np.mean(structure.xyz[0, 1:], axis=0) - structure.xyz[0, 0]
            rmat = rotmat(Vector(*cur_orientation), Vector(*first_bond))
            structure.xyz[0] = (np.matmul(rmat, structure.xyz[0].T)).T
        bond = template.xyz[0, cur_h] - template.xyz[0, cur_n]
        location = template.xyz[0, cur_n] + bond * 1.5
        translate = location - structure.xyz[0, 0]
        structure.xyz[0] += translate
        t, b = structure.topology.to_dataframe()
        t['resSeq'] = i
        structure = md.Trajectory(structure.xyz, md.Topology.from_dataframe(t, b))
        backbone = backbone.stack(structure)
        table, bonds = backbone.topology.to_dataframe()
        resname = list(table.loc[table['chainID'] == 1, 'resName'])[0]
        table['chainID'] = 0
        table.loc[table['resSeq'] > nres, 'resSeq'] = i
        table.loc[table['resSeq'] == i, 'resName'] = resname
        table.loc[:, 'serial'] = np.arange(len(table)) + 1
        top = md.Topology.from_dataframe(table, bonds)
        t2, b2 = top.to_dataframe()
        proper_indices = t2['serial'].to_numpy() - 1
        backbone = md.Trajectory(backbone.xyz[0, proper_indices], top)
        top = backbone.topology
        bond_n = top.select("resid " + str(i-1) + " and name N")[0]
        bond_r = top.select("resname " + resname + " and not name C and not name N and not name CA and not name O")[0]
        atoms = [atom for atom in top.atoms]
        top.add_bond(atoms[np.where(proper_indices == bond_n)[0][0]], atoms[np.where(proper_indices == bond_r)[0][0]])
        backbone.topology = top        

    i += 1

backbone.save_pdb("combined.pdb")

with open("combined.pdb") as f:
    lines = f.readlines()

with open("output.pdb", "w") as f:
    for line in lines:
        if 'UNK' in line or 'CONECT' in line:
            continue
#         if line.startswith("ATOM"):
#             if 'PRO' in line or 'GLU' in line or 'ASP' in line or 'CYS' in line or 'MET' in line:
#                 l = line[:20] + "d    " + line[25:]
#             else:
#                 l = line[:20] + "p    " + line[25:]
#         else:
#             l = line
        f.write(line)
       

def perturb_angle(molecule, residue_num, sigma):
    sidechain = sidechain_indices[residue_num]
    if len(sidechain) == 0:
        return
    cur_energy = energy_function(calculate_min_dist(molecule, residue_num))
    backbone_n = molecule.topology.select("resid " + str(residue_num) + " and name N")[0]
    n_coord = molecule.xyz[0, backbone_n]
    xyz_copy = np.copy(molecule.xyz)
    sidechain_xyz = np.copy(molecule.xyz[0, sidechain])
    ex_first_atom = np.copy(sidechain_xyz[0])
    sidechain_xyz -= n_coord

    first_atom = sidechain_xyz[0]
    perturbation = np.random.normal(0, sigma, 3)
    new_vector = first_atom + perturbation
    normalized = new_vector * (np.linalg.norm(first_atom) / np.linalg.norm(new_vector))
    rmat = rotmat(Vector(*first_atom), Vector(*normalized))
    sidechain_xyz = (np.matmul(rmat, sidechain_xyz.T)).T
    sidechain_xyz += n_coord
    cur_first_atom = np.copy(sidechain_xyz[0])
    sidechain_xyz += (ex_first_atom - cur_first_atom)
    molecule.xyz[0, sidechain] = sidechain_xyz
    new_energy = energy_function(calculate_min_dist(molecule, residue_num))
    if new_energy > cur_energy:
        prob = np.exp(-1 * kT * (new_energy - cur_energy))
        if np.random.uniform() > prob:
            molecule.xyz = xyz_copy
    
def calculate_min_dist(molecule, residue_num):
    res_list = []
    residue_xyz = molecule.xyz[0, sidechain_indices[residue_num]]
    for i, res in enumerate(molecule.top.residues):
        if 'GLY' not in res.name and i != residue_num:
            res_list.append(i)
    minvals = np.full(len(res_list), np.Inf)
    for i, num in enumerate(res_list):
        coords = np.concatenate((backbone_indices[num], sidechain_indices[num]), dtype=np.int64)
        other_xyz = molecule.xyz[0, coords]
        for coord1 in residue_xyz:
            for coord2 in other_xyz:
                dist = np.linalg.norm(coord2 - coord1)
                if minvals[i] > dist:
                    minvals[i] = dist
    return minvals

def energy_function(mindists):
    return np.sum([x ** (-2) for x in mindists])

sigma = 0.01
molecule = md.load('output.pdb')
for k in range(100):
    for j, letter in enumerate(args.seq):
        if k < num_iters[letter]:
            perturb_angle(molecule, j, sigma)
       

molecule.save_pdb("molecule.pdb")

command = "obminimize molecule.pdb"
subprocess.run(command.split(), stdout=subprocess.PIPE)

with open("molecule.pdb") as f:
    lines = f.readlines()

with open("molecule.pdb", "w") as f:
    for line in lines:
        l = line
        if line.startswith("ATOM"):
            if 'PRO' in line or 'GLU' in line or 'ASP' in line or 'CYS' in line or 'MET' in line:
                str2 = "d    "
            else:
                str2 = "p    "
            if line[13:15] == "N ":
                str1 = "NL "
            elif line[13:15] == "C ":
                str1 = "CLP"
            elif line[13:15] == "O ":
                str1 = "OL "
            else:
                str1 = line[13:16]
            l = line[:13] + str1 + line[16:20] + str2 + line[25:]
        f.write(l)
        
command = "rm output.pdb combined.pdb"
subprocess.run(command.split(), stdout=subprocess.PIPE)
