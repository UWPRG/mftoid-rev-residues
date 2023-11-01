#!/bin/bash

# email on start, end, and abortion
#SBATCH --job-name=AAAAA

#SBATCH --output=out.out
#SBATCH --error=error.out
#SBATCH --partition=andrewferguson-gpu
#SBATCH --nodes=1            # SET NUM NODES 
#SBATCH --ntasks-per-node=1  # SETS NUM MPI RANKS (1 PER GPU)
#SBATCH --cpus-per-task=10    # SET NUM THREADS 
#SBATCH --account=pi-andrewferguson
#SBATCH --mem=8GB
#SBATCH --gres=gpu:1
# THIS EXAMPLE USES 1 GPU NODE - 1 MPI TASK - 4 THREADS PER TASK

# SET NUMBER OF MPI TASKS 
# SET NUMBER OF MD STEPS

#LOAD GROMACS MODULE 

NTASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES))

source activate /project2/andrewferguson/berlaga/conda_env
module load gromacs
#module load cuda/10.0
#module load openmpi/2.0.2
#source /project2/andrewferguson/Kirill/plumed_mods/plumed-2.5.2/sourceme.sh
#source /project2/andrewferguson/Kirill/gromacs_2019.2/bin/GMXRC
module load cuda
/project2/andrewferguson/berlaga/bin/gmx_mpi mdrun -ntomp $SLURM_CPUS_PER_TASK --deffnm npt -v -nice 0 -plumed plumed_pb.dat

