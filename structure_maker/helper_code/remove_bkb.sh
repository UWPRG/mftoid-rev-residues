for filename in $(ls noh_residue_pdb/);
do
	NAME1=${filename}
	NAME2="noh_residue_pdb/"
	TOTAL_NAME="$NAME2$NAME1"
       	python nobackbone.py $TOTAL_NAME
	obabel $TOTAL_NAME -O $TOTAL_NAME -d
done

