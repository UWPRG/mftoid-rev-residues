NAME=$1
MINI="NOSOL"
cd structure_maker/
rm molecule*.pdb
python make_structure.py --seq $NAME --mini AD-T --nobackup
cd ../residue_topologies/m_sims
if [ ! -d "$MINI" ]
then 
	mkdir -p "$MINI"
fi 
cp -r ../vacuum_template $MINI/$NAME
cp ../../structure_maker/molecule.pdb $MINI/$NAME/a_insert/
sed -i "s/SEQ/$NAME/g" $MINI/$NAME/mini.sbatch
sed -i "s/AAAAA/$NAME$MINI/g" $MINI/$NAME/4_pb/pbmetad.sbatch
sed -i "s/AAAAA/$NAME$MINI/g" $MINI/$NAME/4_pb/sbatch.sbatch


cd $MINI/$NAME
sbatch --account=pi-andrewferguson mini.sbatch
cd ../../../



