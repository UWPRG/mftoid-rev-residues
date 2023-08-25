for filename in $(ls noh_residue_pdb/);
do
	NAME1=${filename}
	NAME2="noh_residue_pdb/"
	TOTAL_NAME="$NAME2$NAME1"
	final=${NAME1%.pdb}
	st=${NAME1:0:3}
	sed -i "s/$st/$final/g" $TOTAL_NAME
done

