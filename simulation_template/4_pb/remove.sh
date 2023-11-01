rm npt.edr npt.xtc npt.tpr npt.cpt out.out error.out mdout.mdp HILLS* step*
gmx grompp -f mdout_peptoid.mdp -p topol.top -c npt.gro -o npt.tpr -maxwarn 1

