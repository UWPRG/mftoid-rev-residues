import mdtraj as md
import os
import numpy as np  

lines = open('npt.pdb', 'r').readlines()
with open("npt.pdb", "w") as f:
    for line in lines:
        if "SOL" not in line and "CL " not in line and "NA " not in line:
            f.write(line)

