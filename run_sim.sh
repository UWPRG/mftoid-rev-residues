NAME=$1
MINI=$2

cd structure_maker/
rm molecule*.pdb
python make_structure.py --seq $NAME --mini $MINI
cd ../simulations
if [ ! -d "$MINI" ]
then 
	mkdir -p "$MINI"
fi 
cp -r ../simulation_template $MINI/$NAME
cp ../structure_maker/molecule.pdb $MINI/$NAME/a_insert/


cd $MINI/$NAME
python simulation.py
cd ../../../



