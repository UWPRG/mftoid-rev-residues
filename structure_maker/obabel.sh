for filename in $(ls residue_pdb/);
do
	NAME1=${filename}
	NAME2="residue_pdb/"
	TOTAL_NAME="$NAME2$NAME1"
	obabel $TOTAL_NAME -O $TOTAL_NAME 
done

