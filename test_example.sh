NAME=AA
MINI=AD-T

cd structure_maker/
rm molecule*.pdb
python minima.py --seq $NAME --mini $MINI
cd ..
cp -r example_sim test_sim
cp structure_maker/molecule.pdb test_sim/a_insert/
cd test_sim
python simulation.py
cd ../



