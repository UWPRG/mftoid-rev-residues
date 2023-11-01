#!/bin/bash
module unload intelmpi
module load python
#source /project2/andrewferguson/Kirill/plumed_mods/plumed-2.5.2/sourceme.sh
#source /project2/andrewferguson/berlaga/conda_env/bin/GMXRC
module load openbabel
shopt -s expand_aliases

NTASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES))

# SET NUMBER OF OPENMP THREADS
OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

NMOL=1
DIMS=6.0

NTASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES))

NAME=MOL
NAME2=${PWD##*/}

cd a_insert
alias gmx_gpu='/project2/andrewferguson/berlaga/bin/gmx_mpi'
gmx_gpu insert-molecules -ci molecule.pdb -o box.gro -nmol $NMOL -box $DIMS $DIMS $DIMS
gmx_gpu pdb2gmx -f box.gro -p topol.top -o box_mol.gro -water spce -ff charmm36-feb2021

cd ..
cp a_insert/box_mol.gro b_em
cp a_insert/topol.top b_em

cd b_em

sed -i "s/MOLNAME/$NAME/g" npt_peptoid.mdp
sed -i "s/MOLNAME/$NAME/g" nvt_peptoid.mdp

gmx_gpu grompp -f em_peptoid.mdp -c box_mol.gro -p topol.top -o em.tpr -maxwarn 3

gmx_gpu mdrun -ntomp $OMP_NUM_THREADS -v -deffnm em


gmx_gpu grompp -f nvt_peptoid.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 3

gmx_gpu mdrun -ntomp $OMP_NUM_THREADS -v -deffnm nvt

gmx_gpu grompp -f npt_peptoid.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 3

gmx_gpu mdrun -ntomp $OMP_NUM_THREADS -v -deffnm npt

cd ..

cp b_em/npt.gro c_pdb2gmx
cp b_em/npt.cpt c_pdb2gmx
cp b_em/topol.top c_pdb2gmx

cd c_pdb2gmx

gmx_gpu pdb2gmx -f npt.gro -o peptoid.gro -p topol.top -water spce -ff charmm36-feb2021 -ter  

cd ..


cp c_pdb2gmx/peptoid.gro 1_solvate
cp c_pdb2gmx/topol.top 1_solvate
cd 1_solvate

gmx_gpu insert-molecules -ci peptoid.gro -o box.gro -nmol $NMOL -box $DIMS $DIMS $DIMS
gmx_gpu solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

cd ..

cp 1_solvate/solv.gro 2_em
cp 1_solvate/topol.top 2_em

cd 2_em

sed -i "s/MOLNAME/$NAME/g" npt_peptoid.mdp
sed -i "s/MOLNAME/$NAME/g" nvt_peptoid.mdp


echo SOL | gmx_gpu genion -s ions.tpr -o solv.gro -p topol.top -neutral
gmx_gpu grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx_gpu grompp -f em_peptoid.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 3

gmx_gpu mdrun -ntomp $OMP_NUM_THREADS -v -deffnm em


gmx_gpu grompp -f nvt_peptoid.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 3

gmx_gpu mdrun -ntomp $OMP_NUM_THREADS -v -deffnm nvt

gmx_gpu grompp -f npt_peptoid.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 3

gmx_gpu mdrun -ntomp $OMP_NUM_THREADS -v -deffnm npt


cd ..

cp 2_em/npt.gro 3_md
cp 2_em/npt.cpt 3_md
cp 2_em/topol.top 3_md

cd ../3_md

sed -i "s/MOLNAME/$NAME/g" md.mdp

gmx_gpu grompp -f mdout_peptoid.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2

gmx_gpu mdrun -ntomp $OMP_NUM_THREADS -v -deffnm md
echo 1 | gmx_gpu trjconv -f md.xtc -o md_whole.xtc -s md.tpr -pbc whole
echo 1 | gmx_gpu trjconv -f npt.gro -o npt_whole.gro -s md.tpr -pbc whole
gmx_gpu editconf -f npt_whole.gro -o npt_whole.pdb
