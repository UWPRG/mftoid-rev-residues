source /project/andrewferguson/berlaga/conda_env/bin/GMXRC
export OMP_NUM_THREADS=5
module load openmpi/4.1.2+gcc-7.4.0
module load gsl
module load cuda/11.2
module load python

plumed driver --noatoms --plumed plumed_reweight.dat

