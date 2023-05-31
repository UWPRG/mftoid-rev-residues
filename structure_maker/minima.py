import numpy as np
import mdtraj as md
import Bio.PDB
import os
import argparse
import subprocess
import shutil

from Bio.PDB.vectors import Vector, rotaxis, calc_dihedral, rotmat
from Bio.SVDSuperimposer import SVDSuperimposer


# from utils import *
kT = 50

FILENAMES = {
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
MINIMA_FILENAMES = {
    'C7B+T': 'c7beta_plus_trans',
    'C7B-T': 'c7beta_minus_trans',
    'C7B-C': 'c7beta_minus_cis',
    'C7B+C': 'c7beta_plus_cis',
    'AD+T': 'alphaD_plus_trans',
    'AD-T': 'alphaD_minus_trans',
    'AD+C': 'alphaD_plus_cis',
    'AD-C': 'alphaD_minus_cis',
    'A+T': 'alpha_plus_trans',
    'A-T': 'alpha_minus_trans',
    'A+C': 'alpha_plus_cis',
    'A-C': 'alpha_minus_cis'
}

NUM_ITERS = {
    'A': 100,
    'R': 100,
    'N': 100,
    'D': 100,
    'C': 100,
    'E': 100,
    'Q': 100,
    'G': 0,
    'H': 100,
    'J': 100,
    'B': 100,
    'I': 300,
    'L': 100,
    'K': 100,
    'M': 100,
    'F': 300,
    'P': 0,
    'S': 100,
    'T': 300,
    'W': 300,
    'Y': 300,
    'V': 300
}

NTERM_POSITION = np.array([0.0264, -2.7761, 2.590])

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
    "--mini",
    type=str,
    default=None,
    metavar="mini",
    help="minimum configuration",
)
args = parser.parse_args()

if args.seq is None:
    raise ValueError("Sequence Entry Needed")
    
if args.mini is None:
    raise ValueError("Minimum selection needed")



if args.mini not in MINIMA_FILENAMES.keys():
    err_string = str(args.mini) + " is not a valid minimum code. The valid codes are: " + str(minima_filenames.keys())
    raise ValueError(err_string)

mini_string = MINIMA_FILENAMES[args.mini]
backbone = md.load("minima_pdb/" + mini_string + ".pdb")
template = md.load("minima_pdb/" + mini_string + "_nh.pdb")

full_structure = md.Trajectory(np.zeros((1, 0, 3)), md.Topology())
nres = backbone.topology.n_residues
# first_indices = backbone.topology.select("resid 0 or resid 1")
# last_indices = backbone.topology.select("resid 0 or resid 1")
# first_mean = np.mean(backbone.xyz[0, first_indices], axis=0)
# last_mean = np.mean(backbone.xyz[0, last_indices], axis=0)
# linear_chain = last_mean - first_mean
backbone_indices = []
sidechain_indices = []
cur_index = 0
i = 1
for letter in args.seq:
    backbone_indices.append(np.arange(cur_index, cur_index + 4))
    if letter not in FILENAMES.keys():
        err_string = str(letter) + " is not a valid amino acid code."
        raise ValueError(err_string)
    # Read the structure PDB file
    filename = "noh_residue_pdb/" + FILENAMES[letter] + ".pdb"
    cur_index += 4
    if letter == "P":
        # I hate proline!
        orientation = md.load('PROd.pdb', standard_names=False)
        template_backbone_indices = template.topology.select('resSeq ' + str(i) + " and backbone")
        pro_backbone_indices = [3, 1, 0, 2]
        pro_sidechain_indices = [4, 5, 6, 7, 8, 10, 11, 12, 13]
        sup = SVDSuperimposer()
        sup.set(template.xyz[0, template_backbone_indices], orientation.xyz[0, pro_backbone_indices])
        sup.run()
        rms = sup.get_rms()
        rot, tran = sup.get_rotran()
        orientation.xyz[0, pro_sidechain_indices] = np.dot(orientation.xyz[0, pro_sidechain_indices], rot) + tran
#         orientation.superpose(template, atom_indices=pro_backbone_indices, ref_atom_indices=template_backbone_indices)
        structure = md.load(filename, standard_names=False)
        cur_index += structure.topology.n_atoms
        sidechain_indices.append(np.arange(structure.topology.n_atoms) + cur_index)
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
#     elif letter == "G":
#         sidechain_indices.append(np.array([], dtype=np.int32))
#         table, bonds = backbone.topology.to_dataframe()
#         table.loc[table['resSeq'] == i, 'resName'] = "GLY"
#         backbone.topology = md.Topology.from_dataframe(table, bonds)
    else:
        structure = md.load(filename, standard_names=False)
        sidechain_indices.append(np.arange(structure.topology.n_atoms) + cur_index)
        cur_index += structure.topology.n_atoms
        cur_h = template.topology.select("resSeq " + str(i) + " and name H")[0]
        cur_n = template.topology.select("resSeq " + str(i) + " and name N")[0]
        if letter != "A" and letter != "G":
            first_bond = template.xyz[0, cur_h] - template.xyz[0, cur_n]
            cur_orientation = np.mean(structure.xyz[0, 1:], axis=0) - structure.xyz[0, 0]
            rmat = rotmat(Vector(*cur_orientation), Vector(*first_bond))
            structure.xyz[0] = (np.matmul(rmat, structure.xyz[0].T)).T
        bond = template.xyz[0, cur_h] - template.xyz[0, cur_n]
        if letter == 'G':
            location = template.xyz[0, cur_n] + bond
        else:
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
#     cur_first_atom = np.copy(sidechain_xyz[0])
#     sidechain_xyz += (ex_first_atom - cur_first_atom)
    molecule.xyz[0, sidechain] = sidechain_xyz
    new_energy = energy_function(calculate_min_dist(molecule, residue_num))
    if new_energy > cur_energy:
        prob = np.exp(-1 * kT * (new_energy - cur_energy))
        if np.random.uniform() > prob:
            molecule.xyz = xyz_copy
    
def calculate_min_dist(molecule, residue_num):
    res_list = []
    residue_xyz = molecule.xyz[0, sidechain_indices[residue_num]]
    residue_n = molecule.xyz[0, backbone_indices[residue_num][0]]
    residue_ca = molecule.xyz[0, backbone_indices[residue_num][1]]

    for i, res in enumerate(molecule.top.residues):
        if 'GLY' not in res.name and i != residue_num:
            res_list.append(i)
    minvals = np.full(len(res_list), np.Inf)
    for i, num in enumerate(res_list):
        coords = np.concatenate((backbone_indices[num], sidechain_indices[num]), dtype=np.int64)
        other_xyz = molecule.xyz[0, coords]
        for coord1 in residue_xyz:
            if np.any(residue_xyz[0] != coord1) and np.linalg.norm(coord1 - residue_n) < 0.2:
                minvals[i] = 0.02
            elif np.linalg.norm(coord1 - residue_ca) < 0.2:
                minvals[i] = 0.02
            else:
                for coord2 in other_xyz:
                    dist = np.linalg.norm(coord2 - coord1)
                    if minvals[i] > dist:
                        minvals[i] = dist
    return minvals

def energy_function(mindists):
    return np.sum([x ** (-2) for x in mindists])

sigma = 0.01

def minimize_energy():
    molecule = md.load('output.pdb', standard_names=False)
    for k in range(900):
        for j, letter in enumerate(args.seq[:-1]):
            if args.mini == 'A+T' or args.mini == 'A-T':
                if k < 3 * NUM_ITERS[letter]:
                    perturb_angle(molecule, j, sigma)
            else:
                if k < NUM_ITERS[letter]:
                    perturb_angle(molecule, j, sigma)
    molecule.save_pdb('output.pdb')

#ADD HYDROGENS
def add_hydrogens(file1, file2):
    molecule = md.load(file1, standard_names=False)
    new_mtop = molecule.topology
    for j, letter in enumerate(args.seq):
        sup = SVDSuperimposer()
        mtop = molecule.topology
        residue = [res for res in mtop.residues][j]
        traj_h = md.load("residue_pdb/" + FILENAMES[letter] + ".pdb", standard_names=False)
        htop = traj_h.topology
        atoms = np.array([atom for atom in htop.atoms])
        bonds = [bond for bond in htop.bonds]
        if letter != 'G' and letter != 'P':
            fixed_indices = np.concatenate((new_mtop.select("resid " + str(j) + " and not backbone"), new_mtop.select("resid " + str(j) + " and name N"))).astype(np.int64)
            fixed_coords = molecule.xyz[0, fixed_indices] 
            extra_backbone_indices = htop.select("name NL")
            moving_sidechain_indices = htop.select("not name CLP and not name OL and not name NL and not name CA and not element H")
            moving_coords = traj_h.xyz[0, np.concatenate((moving_sidechain_indices, extra_backbone_indices))]
            hydrogen_indices = htop.select("element H and not name HA1 and not name HA2")
            bonded_atom_coords = []
            hydrogens = atoms[hydrogen_indices]
            for atom in hydrogens:
                for bond in bonds:
                    if bond[1].index == atom.index:                 
                        bonded_atom_coords.append(molecule.xyz[0, new_mtop.select("resid " + str(j) + " and name " + bond[0].name)[0]])
                        break 
            bonded_atom_coords = np.array(bonded_atom_coords)
            sup.set(fixed_coords, moving_coords)
            sup.run()
            rms = sup.get_rms()
            rot, tran = sup.get_rotran()
            hydrogen_coords = np.dot(traj_h.xyz[0, hydrogen_indices], rot) + tran
            bonds = hydrogen_coords - bonded_atom_coords
            bls = np.linalg.norm(bonds, axis=1)
            hydrogen_coords = bonded_atom_coords + bonds * (np.divide(0.11, bls)[:, np.newaxis])
            hydrogen_names = [atom.name for atom in hydrogens]
            for name in hydrogen_names:
                mtop.add_atom(name, md.element.hydrogen, residue)
            table, bonds = mtop.to_dataframe()
            new_mtop = md.Topology.from_dataframe(table, bonds)
            new_h_indices = new_mtop.select("resSeq " + str(j+1) + " and element H")
            temp = molecule.xyz[0]
    #         print(new_h_indices)
    #         print(hydrogen_coords)
            for num, idx in enumerate(new_h_indices):
                temp = np.insert(temp, idx, hydrogen_coords[num], axis=0)
            molecule.xyz = np.array([temp])
        fixed_backbone_indices = new_mtop.select("resid " + str(j) + " and backbone and not name O")
        fixed_backbone = molecule.xyz[0, fixed_backbone_indices]
        moving_backbone_indices = htop.select("name NL or name CA or name CLP")
#         print(moving_backbone_indices)
        
        
        moving_backbone = traj_h.xyz[0, moving_backbone_indices]
        ha = htop.select("name HA1 or name HA3 or name HA2 or name HA")
        sup = SVDSuperimposer()
        sup.set(fixed_backbone, moving_backbone)
        sup.run()
        rms = sup.get_rms()
        rot, tran = sup.get_rotran()
        ha_coords = np.dot(traj_h.xyz[0, ha], rot) + tran
        ca_coord = fixed_backbone[1]
        for h in ha_coords:
            bond = ca_coord - h
            bl = np.linalg.norm(bond)
            h = ca_coord + bond * (0.11 / bl)
        temp = molecule.xyz[0]
        if letter == "P":
            mtop.add_atom("HA", md.element.hydrogen, residue)
            table, bonds = mtop.to_dataframe()
            new_mtop = md.Topology.from_dataframe(table, bonds)
            new_ha_index = new_mtop.select("resSeq " + str(j+1) + " and name HA")
            if new_ha_index >= len(molecule.xyz[0]):
                temp = np.append(temp, ha_coords, axis=0)
            else:
                temp = np.insert(temp, new_ha_index, ha_coords[0], axis=0)
        else:
            mtop.add_atom("HA1", md.element.hydrogen, residue)
            mtop.add_atom("HA2", md.element.hydrogen, residue)
            table, bonds = mtop.to_dataframe()
            new_mtop = md.Topology.from_dataframe(table, bonds)
            new_ha1_index = new_mtop.select("resSeq " + str(j+1) + " and name HA1")
            new_ha2_index = new_mtop.select("resSeq " + str(j+1) + " and name HA2")
            if new_ha1_index >= len(molecule.xyz[0]):
                temp = np.append(temp, ha_coords, axis=0)
            else:
                temp = np.insert(temp, new_ha1_index, ha_coords[0], axis=0)
                temp = np.insert(temp, new_ha2_index, ha_coords[1], axis=0)
        molecule.xyz = np.array([temp])
    molecule.save_pdb(file2)

def attach_ace_nme():
    molecule = md.load("molecule.pdb", standard_names=False)
    ace = md.load('residue_pdb/ACEp_f.pdb')
    nme = md.load('residue_pdb/NMEp_f.pdb')
    ace_index = ace.topology.select('name CLP')[0]
    ace_attach = ace.xyz[0, ace_index]
    mol_nterm = NTERM_POSITION
    n_translate = mol_nterm - ace_attach
    ace.xyz[0] += n_translate
    
    n_bond = mol_nterm - molecule.xyz[0, molecule.topology.select("resid 0 and name N")[0]]
    ace.xyz[0] += n_bond * .3
    nme_attach = nme.xyz[0, nme.topology.select('name NL')[0]]
    nres = molecule.topology.n_residues
    template_h = md.load("minima_pdb/" + mini_string + "_h.pdb", standard_names=False)
    if nres == 19:
        mol_cterm = template_h.xyz[0, 134]
    else:
        cterm_idx = template_h.topology.select("resid " + str(nres) + " and name N")[0]
        mol_cterm = template_h.xyz[0, cterm_idx]
    c_translate = mol_cterm - nme_attach
    nme.xyz[0] += c_translate
    c_bond = mol_cterm - molecule.xyz[0, molecule.topology.select("resid " + str(nres - 1) + " and name C")[0]]
    i = nres + 1

    for cap, idx, bond in zip((ace, nme), (0, nres), (n_bond, c_bond)):
        first = np.copy(cap.xyz[0, 0])
        cap.xyz[0] -= first
        cur_orientation = np.mean(cap.xyz[0, cap.topology.select("name H or not element H")[1:]], axis=0)
        rmat = rotmat(Vector(*cur_orientation), Vector(*bond))
        cap.xyz[0] = (np.matmul(rmat, cap.xyz[0].T)).T
        cap.xyz[0] += first
        molecule = molecule.stack(cap)

        table, bonds = molecule.topology.to_dataframe()
        table.loc[table['chainID'] == 1, 'resSeq'] = i 
        i += 1
        table['chainID'] = 0

#         molecule.topology = md.Topology.from_dataframe(table, bonds)
#         table, bonds = molecule.topology.to_dataframe()
#         resname = list(table.loc[table['chainID'] == 1, 'resName'])[0]
#         table['chainID'] = 0
#         table.loc[:, 'serial'] = np.arange(len(table)) + 1
        top = md.Topology.from_dataframe(table, bonds)
#         t2, b2 = top.to_dataframe()
#         proper_indices = t2['serial'].to_numpy() - 1
        molecule = md.Trajectory(molecule.xyz, top)
    molecule.save_pdb('molecule.pdb')


minimize_energy()
# command = "obabel output.pdb -O output2.pdb -d --minimize --steps 1000 --sd"
# subprocess.run(command.split(), stdout=subprocess.PIPE)
add_hydrogens("output.pdb", "molecule.pdb")
nres = attach_ace_nme()

# add_hydrogens("output2.pdb", "molecule2.pdb")        
# command = "obabel molecule.pdb -O molecule3.pdb -d --minimize --steps 1000 --sd"
# subprocess.run(command.split(), stdout=subprocess.PIPE)

fs = ("molecule.pdb",)

for ff in fs:
    with open(ff) as f:
        lines = f.readlines()
    atom_number = 1
    with open(ff, "w") as f:
        for line in lines:
            if "MODEL" in line:
                f.write(line)
            if "ACE" in line:
                numlen = len(str(atom_number))
                str0 = " " * (5 - numlen) + str(atom_number) + "  "
                atom_number += 1
                str2 = "p    " + "0   "                
                if line[13:16] == "CH3":
                    str1 = "CAY"
                elif line[13:18] == "H1  A":
                    str1 = "HY1"
                elif line[13:18] == "H2  A":
                    str1 = "HY2"
                elif line[13:18] == "H3  A":
                    str1 = "HY3"
                else:
                    str1 = line[13:16]
                l = line[:6] + str0 + str1 + line[16:20] + str2 + line[29:]
                f.write(l)
            
        for line in lines:
            l = line
            if "ENDMDL" in line:
                f.write(line)
                continue
            if "ACE" in line or 'MODEL' in line or 'CONECT' in line:
                continue
            elif "NME" in line:
                numlen = len(str(atom_number))
                str0 = " " * (5 - numlen) + str(atom_number) + "  "
                atom_number += 1
                str2 = "p    "
                nres = len(args.seq)
                if nres > 10:
                    str3 = str(nres + 1) + "  "
                else:
                    str3 = str(nres + 1) + "   "
                if line[13:19] == "C   NM":
                    str1 = "CAT"
                elif line[13:15] == "H ":
                    str1 = "HNT"
                elif line[13:18] == "H1  N":
                    str1 = "HT1"
                elif line[13:18] == "H2  N":
                    str1 = "HT2"
                elif line[13:18] == "H3  N":
                    str1 = "HT3"
                else:
                    str1 = line[13:16]
                if "TER" in line:
                    l = line[:6] + str0 + str1 + line[16:20] + str2 + str3 + "\n"
                else:
                    l = line[:6] + str0 + str1 + line[16:20] + str2 + str3 + line[29:]

            elif line.startswith("ATOM"):
                num_str = line[25:27]
                cur_resnum = int(num_str)
                cur_residue = " " + FILENAMES[args.seq[cur_resnum - 1]] + "    "
                numlen = len(str(atom_number))
                str0 = " " * (5 - numlen) + str(atom_number) + " "
                atom_number += 1
#                 if 'PRO' in line or 'GLU' in line or 'ASP' in line or 'CYS' in line or 'MET' in line:
#                     str2 = "d    "
#                 else:
#                     str2 = "p    "
                
                if line[13:15] == "N ":
                    str1 = " NL "
                elif line[13:15] == "C ":
                    str1 = " CLP"
                elif line[13:15] == "O ":
                    str1 = " OL "             
                elif line[13:18] == "H   G":
                    str1 = " HN "
                elif line[13:18] == "HA3 G":
                    str1 = " HA1"         
                else:
                    str1 = line[12:16]
                l = line[:6] + str0 + str1 + cur_residue + line[25:]
            f.write(l)
        
# command = "obabel molecule.pdb -O molecule.pdb"
# subprocess.run(command.split(), stdout=subprocess.PIPE)
command = "rm output.pdb combined.pdb"
subprocess.run(command.split(), stdout=subprocess.PIPE)
