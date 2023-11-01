gmx mdrun -ntomp $SLURM_CPUS_PER_TASK --deffnm npt -cpi npt.cpt -nice 0 -plumed plumed_meta.dat
