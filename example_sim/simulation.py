import subprocess 
import os

"""
If you do not use SLURM as your job scheduler and want this code to work in parallel, please reach out! Or you can simply hardcode the variables NTASKS and CPUS_PER_TASK.

Can also edit NMOL and DIMS
"""


CPUS_PER_TASK = os.environ.get("SLURM_CPUS_PER_TASK", 1)
NMOL = 1
DIMS = 6.0

commands = ["bash", "mini_sim.sh", str(CPUS_PER_TASK), str(NMOL), str(DIMS)]
        
subprocess.run(commands)