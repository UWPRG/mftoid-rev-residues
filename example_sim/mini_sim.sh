#!/bin/bash


# SET NUMBER OF OPENMP THREADS
OMP_NUM_THREADS=$1

NMOL=$2
DIMS=$3


cd a_insert
gmx_mpi insert-molecules -ci molecule.pdb -o box.gro -nmol $NMOL -box $DIMS $DIMS $DIMS
gmx_mpi pdb2gmx -f box.gro -p topol.top -o box_mol.gro -water spce -ff charmm36-feb2021 -v

cd ..
cp a_insert/box_mol.gro b_em
cp a_insert/topol.top b_em

cd b_em

gmx_mpi grompp -f em_peptoid.mdp -c box_mol.gro -p topol.top -o em.tpr -maxwarn 3

gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -v -deffnm em

cd ..

cp b_em/em.gro c_pdb2gmx
cp b_em/topol.top c_pdb2gmx

cd c_pdb2gmx

gmx_mpi pdb2gmx -f em.gro -o peptoid.gro -p topol.top -water spce -ff charmm36-feb2021 -ter  

cd ..


cp c_pdb2gmx/peptoid.gro 1_solvate
cp c_pdb2gmx/topol.top 1_solvate
cd 1_solvate

gmx_mpi solvate -cp peptoid.gro -cs spc216.gro -o solv.gro -p topol.top

cd ..

cp 1_solvate/solv.gro 2_em
cp 1_solvate/topol.top 2_em

cd 2_em

gmx_mpi grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo SOL | gmx_mpi genion -s ions.tpr -o solv.gro -p topol.top -neutral
gmx_mpi grompp -f em_peptoid.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 3

gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -v -deffnm em


gmx_mpi grompp -f nvt_peptoid.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 3

gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -v -deffnm nvt

gmx_mpi grompp -f npt_peptoid.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 3

gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -v -deffnm npt


cd ..

cp 2_em/npt.gro 3_md
cp 2_em/npt.cpt 3_md
cp 2_em/topol.top 3_md

cp 2_em/topol.top 4_pb
cp 2_em/npt.gro 4_pb
cp 2_em/npt.cpt 4_pb

cd 3_md

sed -i "s/MOLNAME/$NAME/g" md.mdp

gmx_mpi grompp -f mdout_peptoid.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2

gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -v -deffnm md
echo 1 | gmx_gpu trjconv -f md.xtc -o md_whole.xtc -s md.tpr -pbc whole
echo 1 | gmx_gpu trjconv -f npt.gro -o npt_whole.gro -s md.tpr -pbc whole
gmx_mpi editconf -f npt_whole.gro -o npt_whole.pdb
cd ../4_pb
gmx_mpi editconf -f npt.gro -o npt.pdb
python remove_solvent.py
python gen_plumed.py
gmx_mpi grompp -f mdout_peptoid.mdp -c npt.gro -o npt.tpr -maxwarn 2
gmx_mpi mdrun -ntomp $OMP_NUM_THREADS --deffnm npt -nice 0 -plumed plumed_meta.dat
echo 1 | gmx_mpi trjconv -f npt.xtc -o npt_whole.xtc -s npt.tpr -pbc whole
echo 1 | gmx_mpi trjconv -f npt.gro -o npt_whole.gro -s npt.tpr -pbc whole
gmx_mpi editconf -f npt_whole.gro -o npt_whole.pdb
rm -f colvar_reweight.dat
bash reweight.sh
