import numpy as np
import copy
import mdtraj as md
import os
import argparse
import subprocess
import sys
from openbabel import pybel
from Bio.PDB import PDBIO
from Bio.PDB.vectors import Vector, rotaxis, calc_dihedral, rotmat
from Bio.SVDSuperimposer import SVDSuperimposer
from PeptideBuilder import Geometry
import PeptideBuilder
from _gromacs_utils import *
import shutil
from collections import defaultdict
import itertools
import json

MINIMA_FILENAMES = ['AD-C', 'AD+C', 'A-C', 'A+C', 'C7B-C', 'C7B+C', 'AD-T', 'AD+T', 'A-T', 'A+T', 'C7B-T', 'C7B+T', 'HELIX']
MINIMA_ANGLES =  [[ -82,  175,    0],
       [  82, -176,    0],
       [ -61,  -50,  -11],
       [  63,   49,   11],
       [-135,   74,   -3],
       [ 128,  -71,    2],
       [ -82,  174,  180],
       [  83, -177, -180],
       [ -59,  -64,  176],
       [  59,   64, -177],
       [ -82,   80, -177],
       [ 105, -185,  178],
       [ -60,  -50,  180]]

MINIMA_DICT = dict(zip(MINIMA_FILENAMES, MINIMA_ANGLES))

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
    'V': 'VALp',
    '1': 'NTBp',
    '2': 'NPHp',
    '3': 'BTM',
    '4': 'BTE',
    '5': 'CTM',
    '6': 'CTE',
    '7': 'FTM',
    '8': 'FTE',
    '9': 'ITM',
    '0': 'ITE',
    'Z': 'NSPE',
    'X': 'NAEk',
    '@': 'EME',
    '?': 'NOE',
    'U': 'DECp'
}

endcap_params = {
    "nterm" : 'ACEp_f', 
    "cterm" : 'NMEp_f', 
    "n_indexname" : 'CLP', 
    "c_indexname" : 'NL'
}

NTERM_POSITION = np.zeros(3)
CTERM_POSITION = np.zeros(3)

def check_files(necessary_files):
    for file in necessary_files:
        if not os.path.isfile(file):
            raise FileNotFoundError(file + "does not exist.")

def get_dihedral(traj, dihedral, i):
    dihedral_list = ["phi", "psi", "omega"]
    atom_angles = [["CLP", "NL", "CA", "CLP"], ["NL", "CA", "CLP", "NL"], ["CA", "CLP", "NL", "CA", "CAT"]]
    adjustment = [[-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1, 1]]
    k = dihedral_list.index(dihedral)
    at_angles = atom_angles[k]
    res_nums = np.array(adjustment[k]) + i
    search_strings = ["resid " + str(rn) + " and name " + str(ang) for rn, ang in zip(res_nums, at_angles)]
    search_indices = traj.topology.select(" or ".join(search_strings))
    return md.compute_dihedrals(traj, [search_indices])
    
def make_polyglycine_angles(sequence, angles):
    global NTERM_POSITION
    global CTERM_POSITION

    ##MAKE POLYGLYCINE 
    if angles.ndim == 1:
        assert len(angles) == 3
        all_angles = np.tile(angles, (len(sequence) + 2, 1))
    else:
        assert angles.ndim == 2
        assert angles.shape[0] == len(sequence) + 2
        assert angles.shape[1] == 3
        all_angles = angles

    new_seq = "G" + sequence + "G"
    for i in range(len(sequence) + 2):
        if new_seq[i] == 'P':
            geo = Geometry.geometry("P")
        else:
            geo = Geometry.geometry("G")
        geo.phi = all_angles[i, 0]
        geo.psi_im1 = all_angles[i, 1]
        geo.omega = all_angles[i, 2]

        if i == 0:
            structure = PeptideBuilder.initialize_res(geo)
        else:
            PeptideBuilder.add_residue(structure, geo)

    out = PDBIO()
    out.set_structure(structure)
    out.save("template.pdb")

    #SET CTERM/NTERM POSITIONS

    template = md.load("template.pdb", standard_names=False)
    NTERM_POSITION = template.xyz[0, template.top.select("resSeq 1 and name C")[0]]
    CTERM_POSITION = template.xyz[0, template.top.select("resSeq {} and name N".format(len(sequence) + 2))[0]]


    #MAKE -H AND -NH STRUCTURES
    command = "obabel template.pdb -O template_h.pdb -h"
    subprocess.run(command.split(), stdout=subprocess.PIPE)
    template_h = md.load('template_h.pdb', standard_names=False)
    top = copy.deepcopy(template_h.top)
    carbon_top = copy.deepcopy(template_h.top)
    onswitch_nitrogen = np.zeros(top.n_atoms, dtype=np.bool_)
    onswitch_carbon = np.zeros(top.n_atoms, dtype=np.bool_)

    backbone_idxs = top.select('backbone or element C')
    onswitch_nitrogen[backbone_idxs] = 1
    onswitch_carbon[backbone_idxs] = 1

    nterm_coord_exists = False
    ha_names = ['HA1', 'HA2', '', '']
    ha_counters = np.zeros(len(new_seq), dtype=np.int32)
    pro_counters = [defaultdict(int) for k in range(len(new_seq))]
    for bond in template_h.top.bonds:
        if bond[1].element.name == 'hydrogen':
            if bond[0].element.name == 'nitrogen':
                if bond[0].residue.resSeq == 1 and not nterm_coord_exists:
                    nterm_coord_exists = True
                else:
                    onswitch_nitrogen[bond[1].index] = 1
                    if new_seq[bond[0].residue.index] == 'G':
                        onswitch_carbon[bond[1].index] = 1
                        carbon_top.atom(bond[1].index).name = 'HN'
            if bond[0].name == 'CA':
                onswitch_carbon[bond[1].index] = 1
                onswitch_nitrogen[bond[1].index] = 1
                if bond[0].residue.name == 'PRO':
                    carbon_top.atom(bond[1].index).name = 'HA'
                else:
                    carbon_top.atom(bond[1].index).name = ha_names[ha_counters[bond[0].residue.index]]
                ha_counters[bond[0].residue.index] += 1
            elif bond[0].residue.name == 'PRO' and len(bond[0].name) > 1:
                onswitch_carbon[bond[1].index] = 1
                onswitch_nitrogen[bond[1].index] = 1
                carbon_top.atom(bond[1].index).name = "H" + bond[0].name[1] + str(pro_counters[bond[0].residue.index][bond[0].name] + 1)
                pro_counters[bond[0].residue.index][bond[0].name] += 1
    
    for i, (ha_counter, pro_counter) in enumerate(zip(ha_counters, pro_counters)):
        if new_seq[i] == 'P':
            if ha_counter != 1:
                print(f"Proline should have 1 HA. You have {ha_counter}. OpenBabel could not place HA correctly due to aphysical bond angles. Try a different dihedral configuration.")
                raise ValueError(f"Proline should have 1 HA. You have {ha_counter}. OpenBabel could not place HA correctly due to aphysical bond angles. Try a different dihedral configuration.")
            for key, counter in pro_counter.items():
                print(f"Your peptoid should have 2 {key}s. You have {counter}. OpenBabel could not place {key} correctly due to aphysical bond angles. Try a different dihedral configuration.")
                if counter != 2:
                    raise ValueError(f"Your peptoid should have 2 {key}s. You have {counter}. OpenBabel could not place {key} correctly due to aphysical bond angles. Try a different dihedral configuration.")
        else:
            if ha_counter != 2:
                print(f"Your amino acid {FILENAMES[new_seq[i]]} should have 2 HAs. You have {ha_counter}. OpenBabel could not place HAs correctly due to aphysical bond angles. Try a different dihedral configuration.")
                raise ValueError(f"Your amino acid should have 2 HAs. You have {ha_counter}. OpenBabel could not place HAs correctly due to aphysical bond angles. Try a different dihedral configuration.")
    for i in range(top.n_atoms - 1, -1, -1):
        if not onswitch_nitrogen[i]:
            top.delete_atom_by_index(i)
    for i in range(carbon_top.n_atoms - 1, -1, -1):
        if not onswitch_carbon[i]:
            carbon_top.delete_atom_by_index(i)

    new_traj = md.Trajectory(template_h.xyz[0, onswitch_carbon], carbon_top)
    real = new_traj.atom_slice(new_traj.top.select("resSeq 2 to " + str(len(sequence) + 1)))    
    real.save_pdb("template.pdb")
    
    new_traj = md.Trajectory(template_h.xyz[0, onswitch_nitrogen], top)
    real_nh = new_traj.atom_slice(new_traj.top.select("resSeq 2 to " + str(len(sequence) + 1)))
    real_nh.save_pdb("template_nh.pdb")
    
    
    for suffix in ["", "_h", "_nh"]:
        for resid in range(2, len(sequence) + 2):
            aa_str = 'GLY'
            if new_seq[resid - 1] == 'P':
                aa_str = 'PRO'
            if resid > 10:
                commands = 'sed|-i|s/{} A  {}/{} A  {}/g|template{}.pdb'.format(aa_str, resid, aa_str, resid - 1, suffix)
            elif resid == 10:
                commands = 'sed|-i|s/{} A  10/{} A   9/g|template{}.pdb'.format(aa_str, aa_str, suffix)
            else:
                commands = 'sed|-i|s/{} A   {}/{} A   {}/g|template{}.pdb'.format(aa_str, resid, aa_str, resid - 1, suffix)
            subprocess.run(commands.split('|'), stdout=subprocess.PIPE)

def make_peptide_angles():
    global NTERM_POSITION
    global CTERM_POSITION
    if angles.ndim == 1:
        assert len(angles) == 3
        all_angles = np.tile(angles, (len(sequence) + 2, 1))
    else:
        assert angles.ndim == 2
        assert angles.shape[0] == len(sequence) + 2
        assert angles.shape[1] == 3
        all_angles = angles

    new_seq = "G" + sequence + "G"
    for i in range(len(sequence) + 2):
        geo = Geometry.geometry(new_seq[i])
        geo.phi = all_angles[i, 0]
        geo.psi_im1 = all_angles[i, 1]
        geo.omega = all_angles[i, 2]

        if i == 0:
            structure = PeptideBuilder.initialize_res(geo)
        else:
            PeptideBuilder.add_residue(structure, geo)
    
    out = PDBIO()
    out.set_structure(structure)
    out.save("template.pdb")
    command = "obabel template.pdb -O template_h.pdb -h"
    subprocess.run(command.split(), stdout=subprocess.PIPE)
    template = md.load("template.pdb", standard_names=False)
    template_h = md.load('template_h.pdb', standard_names=False)
    NTERM_POSITION = template.xyz[0, template.top.select("resSeq 1 and name C")[0]]
    CTERM_POSITION = template.xyz[0, template.top.select("resSeq {} and name N".format(len(sequence) + 2))[0]]
    real_h = template_h.atom_slice(template_h.top.select("resSeq 2 to " + str(len(sequence) + 1)))
    real_h.save_pdb("template.pdb")
    os.remove("template_h.pdb")
    return real_h
    
def load_sequence(sequence, filename):
    
    backbone = md.load("template.pdb", standard_names=False)
    template = md.load("template_nh.pdb", standard_names=False)
   
    nres = backbone.topology.n_residues
    i = 1
    for letter in sequence:
        # backbone_indices.append(np.arange(cur_index, cur_index + 4))
        if letter not in FILENAMES.keys():
            err_string = str(letter) + " is not a valid amino acid code."
            raise ValueError(err_string)
        
        # Read the structure PDB file
        filename2 = ".nobkb_residue_pdb/" + FILENAMES[letter] + ".pdb"
        
        if letter != 'P' and letter != 'G':
            structure = md.load(filename2, standard_names=False)
            

            #Place sidechain where the N-hydrogen used to be
            cur_h = template.topology.select("resSeq " + str(i) + " and name H")[0]
            cur_n = template.topology.select("resSeq " + str(i) + " and name N")[0]
            if letter != "A" and letter != "G":
                #Adjust orientation of side chain, pointing away from the backbone to avoid clash
                first_bond = template.xyz[0, cur_h] - template.xyz[0, cur_n]
                cur_orientation = np.mean(structure.xyz[0, 1:], axis=0) - structure.xyz[0, 0]
                rmat = rotmat(Vector(*cur_orientation), Vector(*first_bond))
                structure.xyz[0] = (np.matmul(rmat, structure.xyz[0].T)).T
            bond = template.xyz[0, cur_h] - template.xyz[0, cur_n]
            #increase C--N bond length to a normal C--N bond
            if letter == 'G':
                location = template.xyz[0, cur_n] + bond
            else:
                location = template.xyz[0, cur_n] + bond * 1.5
            translate = location - structure.xyz[0, 0]
            structure.xyz[0] += translate

            #Adjust atom indices via MDTraj's dataframe capabilities

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
            
        i += 1
    
    backbone.save_pdb(filename)

    with open(filename, "r") as f:
        lines = f.readlines()

    with open(filename, "w") as f:
        #remove unnecessary bonds from file
        for line in lines:
            if 'UNK' in line or 'CONECT' in line:
                continue
            f.write(line)
   

def _osremove(f):
    if os.path.isfile(f):
        os.remove(f)

def minimize_energy(sequence, filename_before, filename_after, angles):
    for file in os.listdir('.gmx_files'):
        _osremove(os.path.join('.gmx_files', file))
    _osremove("posre.itp")    
    _osremove("#posre.itp.1#")
    
    molecule = md.load(filename_before, standard_names=False)
    check_files([filename_before])  
    insert_molecules_cubic(True, pdb=filename_before, box=".gmx_files/box.gro", dim=max(float(len(sequence)), 3.5))
    _osremove(filename_before)
    check_files([".gmx_files/box.gro"])
    pdb2gmx(True, input_pdb=".gmx_files/box.gro", output_gro=".gmx_files/struct.gro", top=".gmx_files/topol.top", ff="charmm36-feb2021", np=1)
    check_files([".gmx_files/struct.gro"])

    editconf(True, ".gmx_files/struct.gro", ".gmx_files/centered.gro", center=True)
    # fix_dihedral_posres(molecule, angles)
    fix_posres(molecule, 0, molecule.top.n_residues - 1)
    check_files([".em_peptoid.mdp", ".gmx_files/centered.gro", ".gmx_files/topol.top", ".gmx_files/posre.itp"])
    grompp(True, mdp=".em_peptoid.mdp", coord=".gmx_files/centered.gro", top=".gmx_files/topol.top", tpr=".gmx_files/em.tpr", use_index=False, maxwarn=3, np=1)
    check_files([".gmx_files/em.tpr"])
    mdrun(True, deffnm=".gmx_files/em", plumed=False, np=1)
    fix_posres(molecule, 0, molecule.top.n_residues, 20000)
    check_files([".nvt_peptoid.mdp", ".gmx_files/em.gro", ".gmx_files/centered.gro", ".gmx_files/topol.top", "posre.itp"])
    grompp(True, mdp=".nvt_peptoid.mdp", coord=".gmx_files/em.gro", top=".gmx_files/topol.top", tpr=".gmx_files/nvt.tpr", use_index=False, maxwarn=3, np=1)
    check_files([".gmx_files/nvt.tpr"])
    mdrun(True, deffnm=".gmx_files/nvt", plumed=False, np=1)
    check_files([".gmx_files/em.gro"])
    editconf(True, ".gmx_files/em.gro", filename_after, center=False)
    _osremove("posre.itp")    
    _osremove("#posre.itp.1#")
    
def fix_posres(traj, res_begin, res_end, const=100000):
    backbone_idxs = traj.topology.select(f"resSeq {res_begin} to {res_end} and (name CA or name NL or name CLP)")
    with open('.gmx_files/posre.itp', "w") as f:
        f.write('[ position_restraints ]\n\n')
        for d in backbone_idxs:
            f.write(f'{d+1}\t1\t{const}\t{const}\t{const}\n')
            

def fix_dihedral_posres(traj, angles):
    dihedral_list = ["phi", "psi", "omega"]
    atom_angles = [["CLP", "NL", "CA", "CLP"], ["NL", "CA", "CLP", "NL"], ["CA", "CLP", "NL", "CA", "CAT"]]
    adjustment = [[-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1, 1]]
    with open('.gmx_files/posre.itp', "w") as f:
        f.write('[ dihedral_restraints ]\n\n')
        for k, d in enumerate(dihedral_list):
            k = dihedral_list.index(d)
            at_angles = atom_angles[k]
            for i in range(1, traj.top.n_residues - 1):
                res_nums = np.array(adjustment[k]) + i
                if not (np.any(res_nums < 0) or np.any(res_nums > traj.top.n_residues)):
                    search_strings = ["resid " + str(rn) + " and name " + str(ang) for rn, ang in zip(res_nums, at_angles)]
                    search_indices = traj.topology.select(" or ".join(search_strings)) + 1
                    if len(search_indices) == 4:
                        lst = [str(s) for s in search_indices]
                        lst.extend(["9", str(angles[k]), "0", "80"])
                        f.write("\t".join(lst) + '\n')


def attach_endcap(filename, nterm, cterm, n_indexname, c_indexname):
    #Attach end caps to molecule 
    molecule = md.load(filename, standard_names=False)
    ace = md.load(f'.nobkb_residue_pdb/{nterm}.pdb', standard_names=False)
    nme = md.load(f'.nobkb_residue_pdb/{cterm}.pdb', standard_names=False)
    ace_index = ace.topology.select(f'name {n_indexname}')[0]
    ace_attach = ace.xyz[0, ace_index]
    n_translate =  NTERM_POSITION - ace_attach
    ace.xyz[0] += n_translate    
    n_bond = NTERM_POSITION - molecule.xyz[0, molecule.topology.select("resid 0 and name N")[0]]
    nme_attach = nme.xyz[0, nme.topology.select(f'name {c_indexname}')[0]]
    nres = molecule.topology.n_residues
    c_translate = CTERM_POSITION - nme_attach
    nme.xyz[0] += c_translate
    c_bond = CTERM_POSITION - molecule.xyz[0, molecule.topology.select("resid " + str(nres - 1) + " and name C")[0]]
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

        top = md.Topology.from_dataframe(table, bonds)

        molecule = md.Trajectory(molecule.xyz, top)
    molecule.save_pdb(filename)

def pdb_to_peptoid_forcefield(sequence, filename):
    with open(filename) as f:
        lines = f.readlines()
    atom_number = 1

    with open(filename, "w") as f:
        #Adjust PDB file so it includes peptoid atom names, rather than peptide atom names for the caps and backbone. 
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
            if "ACE" in line or 'MODEL' in line or 'CONECT' in line or 'CRYST' in line or 'REMARK' in line:
                continue
            elif "NME" in line:
                numlen = len(str(atom_number))
                str0 = " " * (5 - numlen) + str(atom_number) + "  "
                atom_number += 1
                str2 = "p   "
                nres = len(sequence)
                if nres > 10:
                    str3 = str(nres + 1) + "   "
                else:
                    str3 = " " + str(nres + 1) + "   "
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
                num_str = line[24:26]
                cur_resnum = int(num_str)
                cur_residue = " " + FILENAMES[sequence[cur_resnum - 1]] + "   "
                numlen = len(str(atom_number))
                str0 = " " * (5 - numlen) + str(atom_number) + " "
                atom_number += 1

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
                l = line[:6] + str0 + str1 + cur_residue + line[24:]
            f.write(l)




def create_peptoid(sequence, angle_seq, filename="molecule.pdb", endcap_par=endcap_params):
    if not filename.endswith('.pdb'):
        filename = filename + '.pdb'
    make_polyglycine_angles(sequence, angle_seq)
    load_sequence(sequence, "molecule.pdb")
    attach_endcap("molecule.pdb", **endcap_par)
    pdb_to_peptoid_forcefield(sequence, "molecule.pdb")
    minimize_energy(sequence, "molecule.pdb", filename, angle_seq)
    
    _osremove('template.pdb')
    _osremove('template_h.pdb')
    _osremove('template_nh.pdb')
    
    return md.load(filename, standard_names=False)

def create_peptide(sequence, angle_seq, filename="peptide.pdb", endcap_par=endcap_params):
    make_peptide_angles(sequence, angle_seq)
    attach_endcap(filename, **endcap_par)
    minimize_energy(sequence, "molecule.pdb", filename, angle_seq)
    
def _run_tide_or_toid(tide_or_toid):
    if tide_or_toid:
        return create_peptide
    else:
        return create_peptoid

def run():
    parser = argparse.ArgumentParser(
        description="Molecule Inserter for peptoids"
    )
    parser.add_argument(
        "--seq",
        type=str,
        default=None,
        metavar="seq",
        help="Input Sequence"
    )
    parser.add_argument(
        "--mini",
        type=str,
        default=None,
        metavar="mini",
        help="minimum configuration"
    )
    parser.add_argument(
        "--phi",
        type=int,
        default=None,
        metavar="phi",
        help="minimum configuration"
    )
    parser.add_argument(
        "--psi",
        type=int,
        default=None,
        metavar="psi",
        help="minimum configuration"
    )
    parser.add_argument(
        "--omega",
        type=int,
        default=None,
        metavar="omega",
        help="minimum configuration"
    )
    parser.add_argument(
        "--anglefile",
        type=str,
        default=None,
        metavar="anglefile",
        help="minimum configuration"
    )
    parser.add_argument(
        "--file",
        type=str,
        default="molecule",
        metavar="file",
        help="filename under which to save final structure"
    )
    parser.add_argument(
        "--tide",
        action="store_true",
        help="Make a peptide, not a peptoid"
    )
    parser.add_argument(
        "--nobackup",
        action="store_true",
        help="Overwrite existing files."
    )
    args = parser.parse_args()

    if args.seq is None:
        raise ValueError("Sequence Entry Needed")

    # Load the minimum that the peptoid is in

    if args.file.endswith(".pdb"):
        fname = args.file
    else:
        fname = args.file + ".pdb"
    path = os.path.join(os.getcwd(), fname)
    num_backups = 0
    if not args.nobackup:
        while os.path.exists(path):
            if num_backups == 0:
                print("WARNING: " + fname[:-4]  + ".pdb already exists. Backing up file to " + fname[:-4] + "1.pdb")
            else:
                print("WARNING: " + fname[:-4] + str(num_backups) + ".pdb already exists. Backing up file to " + fname[:-4] + str(num_backups + 1) + ".pdb")
            num_backups += 1
            path = os.path.join(os.getcwd(), fname[:-4] + str(num_backups) + ".pdb")
        while num_backups > 0:
            if num_backups == 1:
                  command = "cp " + fname + " " + fname[:-4] + "1.pdb"
            else:
                  command = "cp " + fname[:-4] + str(num_backups - 1) + ".pdb" + " " + fname[:-4] + str(num_backups) + ".pdb"
            subprocess.run(command.split(), stdout=subprocess.PIPE)
            num_backups -= 1

    # assert args.mini is not None or (args.phi is not None and args.psi is not None and args.omega is not None) or args.anglefile is not None

    if args.mini is not None:

        # try:
        angles = np.array(MINIMA_DICT[args.mini])
        _run_tide_or_toid(args.tide)(args.seq, angles, fname)
        return
        # except:
        #     print("Ramachandran minimum error -- minimum mist be one of: " + str(MINIMA_FILENAMES) + ". Your minimum is " + args.mini)
    if args.anglefile is not None:
        try:
            angles = np.loadtxt(args.anglefile)
            _run_tide_or_toid(args.tide)(args.seq, angles, fname)
            return
        except:
            print("Angle file error -- Could not load text file as numpy: ", args.anglefile)
    if args.phi is not None and args.psi is not None and args.omega is not None:
        angles = np.array([args.phi, args.psi, args.omega])
        _run_tide_or_toid(args.tide)(args.seq, angles, fname)
        return
    else:
        print("Input Error -- Must enter: minimum (--mini) OR phi/psi/omega triplet (--phi, --psi, --omega) OR text file containing angles (--anglefile)")
        return
run()


            
