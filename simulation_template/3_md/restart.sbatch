#!/bin/bash

# email on start, end, and abortion
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berlaga@uchicago.edu

#SBATCH --job-name=NAME

#SBATCH --output=out.out
#SBATCH --partition=broadwl
#SBATCH --nodes=1            # SET NUM NODES 
#SBATCH --ntasks-per-node=1  # SETS NUM MPI RANKS (1 PER GPU)
#SBATCH --cpus-per-task=10    # SET NUM THREADS 
#SBATCH --account=pi-andrewferguson


# SET NUMBER OF MPI TASKS 
# SET NUMBER OF MD STEPS

#LOAD GROMACS MODULE 


module load cuda/10.0
module load openmpi/2.0.2
source /project2/andrewferguson/Kirill/plumed_mods/plumed-2.5.2/sourceme.sh
source /project2/andrewferguson/Kirill/gromacs_2019.2/bin/GMXRC
module load openbabel

NTASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES))

# SET NUMBER OF OPENMP THREADS
OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpirun -np $NTASKS gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -v -deffnm md -cpi md.cpt
echo 1 | gmx_mpi trjconv -f md.xtc -o md_whole.xtc -s md.tpr -pbc whole
echo 1 | gmx_mpi trjconv -f npt.gro -o npt_whole.gro -s md.tpr -pbc whole
obabel npt_whole.gro -O npt_whole.pdb
