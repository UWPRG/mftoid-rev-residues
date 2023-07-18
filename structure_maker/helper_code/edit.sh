for filename in $(ls minima_pdb/*s.pdb)
do
	sed -i "s/GLY/UNK/g" $filename
done
