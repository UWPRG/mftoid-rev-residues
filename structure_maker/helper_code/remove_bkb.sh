for filename in $(ls residue_pdb/);
do
	NAME1=${filename}
	NAME2="noh_residue_pdb/"
	TOTAL_NAME="$NAME2$NAME1"
       	python nobackbone.py $TOTAL_NAME
	obabel $TOTAL_NAME -d -O $TOTAL_NAME	
done

