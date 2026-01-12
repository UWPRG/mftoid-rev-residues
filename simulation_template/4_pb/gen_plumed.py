import mdtraj as md
import os
import numpy as np  

def stringlist(items, sep=','):
    strlist = [str(a) for a in items]
    return sep.join(strlist)

traj = md.load_pdb("npt.pdb")
top = traj.topology

angles = ["phi", "psi", "omega"]
atom_angles = [["CLP", "NL", "CA", "CLP"], ["NL", "CA", "CLP", "NL"], ["CA", "CLP", "NL", "CA", "CAT"]]
adjustment = [[-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1, 1]]

cvs = []
short_cvs = []
alists = []
rg_atoms = []
for k in range(len(angles)):
    angle = angles[k]
    at_angles = atom_angles[k]
    for i in range(1, top.n_residues - 1):
        res_nums = np.array(adjustment[k]) + i
        search_strings = ["resid " + str(rn) + " and name " + str(ang) for rn, ang in zip(res_nums, at_angles)]
        search_indices = top.select(" or ".join(search_strings))
        for s in search_indices:
            if 1 + s not in rg_atoms:
                rg_atoms.append(1 + s)
        atomlist = stringlist(1 + np.array(search_indices))
        short_cvs.append(angle + str(i))
        if k != 2:
            cvs.append(angle + str(i) + ": TORSION ATOMS=" + atomlist + "\n")
            alists.append(atomlist)
        if k == 2:
            atomlist = stringlist(1 + np.array(search_indices))
            cvs.append(angle + str(i) + ": TORSION ATOMS=" + atomlist + "\n")
            alists.append(atomlist)
rg_atoms = sorted(rg_atoms)                                                              
labels = ['alpha_d_minus',
          'alpha_d_minus_trans',
          'alpha_d_plus', 
          'alpha_d_plus_trans',
          'alpha_plus', 
          'alpha_plus_trans',
          'alpha_minus', 
          'alpha_minus_trans',
          'c7beta_plus', 
          'c7beta_plus_trans',
          'c7beta_minus', 
          'c7beta_minus_trans']
angle_vals = [[-1.57, 3.14, 0],
              [-1.57, 3.14, 3.14],
              [1.57, 3.14, 0],
              [1.57, 3.14, 3.14],
              [0.87, 1.00, 0],
              [0.87, 1.00, 3.14],
              [-0.87, 1.00, 0],
              [-0.87, 1.00, 3.14],
              [2.09, -1.30, 0],
              [2.09, -1.30, 3.14],
              [-2.09, 1.30, 0],
              [-2.09, 1.30, 3.14]]

n_triplets = len(alists) / 3

with open("plumed_meta.dat", "w") as f:
    # beginning line
    f.write("#! vim:ft=plumed\n")
    f.write("rg: GYRATION TYPE=RADIUS ATOMS=" + stringlist(rg_atoms) + "\n")
    # alphabetas
    for label, agls in zip(labels, angle_vals):
        f.write("# " + label + "\n")
        f.write("ALPHABETA ...\n")
        for a_idx in range(len(alists)):
            f.write("ATOMS" + str(a_idx + 1) + "=" + alists[a_idx] + " REFERENCE" + str(a_idx + 1) + "=" + str(agls[int(a_idx // n_triplets)]) + "\n")
        f.write("LABEL=" + label + "\n")
        f.write("... ALPHABETA\n")
    # torsion angles
    f.write("# each dihedral angle for a specific residue#\n")
    for cv in cvs:
        f.write(cv)
    #pbmetad command
    f.write("PBMETAD ...\n")
    scv_rng = ",".join(short_cvs)
    f.write("ARG=rg," + ",".join(labels) + "," + scv_rng + "\n")
    sig_ab_rng = ",".join("1" for n in range(12))
    min_ab_rng = ",".join("-1" for n in range(12))
    max_ab_rng = ",".join("50" for n in range(12))

    sig_rng = ",".join("0.35" for n in range(len(cvs)))
    min_rng = ",".join("-pi" for n in range(len(cvs)))
    max_rng = ",".join("pi" for n in range(len(cvs)))
 
    f.write("SIGMA=0.1," + sig_ab_rng + "," + sig_rng + "\n")
    f.write("HEIGHT=3.0 #kJ/mol\n" + "PACE=500\n" + "BIASFACTOR=20\n" + "LABEL=pb\n")
    f.write("GRID_MIN=0," + min_ab_rng + "," + min_rng + "\n")
    f.write("GRID_MAX=3," + max_ab_rng + "," + max_rng + "\n")
    
    file_ab_rng = ",".join(["HILLS." + l for l in labels])
    file_scv_rng = ",".join(["HILLS." + scv for scv in short_cvs])
    f.write("FILE=HILLS.rg," + file_ab_rng + "," + file_scv_rng + "\n")
    f.write("... PBMETAD\n")
    # print command
    f.write("PRINT ARG=rg," + ",".join(labels) + "," + scv_rng + ",pb.bias STRIDE=500 FILE=COLVAR\n")

with open("plumed_reweight.dat", "w") as f:
    f.write("RESTART\n")
    for scv in short_cvs:
        f.write(scv + ": READ FILE=COLVAR VALUES=" + scv + " IGNORE_TIME IGNORE_FORCES EVERY=1\n")
    f.write("\n")
    f.write("PBMETAD ...\n")
    f.write("ARG=" + scv_rng + "\n")
    f.write("SIGMA=" + sig_rng + "\n")
    f.write("HEIGHT=0 #kJ/mol\n" + "PACE=500000\n" + "TEMP=300.0\n" + "BIASFACTOR=20\n")
    f.write("GRID_MIN="+ min_rng + "\n")
    f.write("GRID_MAX="+ max_rng + "\n")
    f.write("LABEL=pb\n" + "... PBMETAD\n")
    f.write("PRINT ARG=" + scv_rng + ",pb.bias STRIDE=1 FILE=colvar_reweight.dat\n") 

                
