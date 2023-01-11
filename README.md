# mftoid-rev-residues

## Pfaendtner + Ferguson Research Team Collaboration

## GOAL: Build library of residues based on the Weiser+Santiso modified MFTOID force field peptoid backbone

### REPOSITORY CONTENTS
	
1. merged.rtp: file containing atoms, bonds, atomtypes, and charges for each peptoid amino acid mimic

2. ffbonded: directory containing 4 files (bondtypes.itp contains updated bond types, angletypes.itp angle types, dihedraltypes.itp proper dihedral types, and improper.itp the improper dihedral types)

3. ffnonbonded.itp: file containing the carbonyl carbon peptoid parameter from Mirijanian's MFTOID

4. residue_pdb: a directory containing PDBs of each of the peptoid versions of the amino acids, orderedand labeled according to updated forcefield residues found in merged.rtp

### INSTRUCTIONS

Download the CHARMM36 forcefield from the MacKerell Lab website; there's been some updates between February 2021 and July 2022 on forcefield file structure. We have been using the February 2021 version of the FF and older.

To update the downloaded version of your CHARMM FF for peptoid amino acid mimic capabilities, start by copying the information found in this Github repositroy in the merged.rtp file directly into the merged.rtp doc in your CHARMM36 FF file.

Next, add the peptoid parameter TC found here in ffnonbonded.itp into the ffnonbonded.itp doc in your CHARMM36 FF file.

Next, add the information contained here in bondtypes.itp, angletypes.itp, dihedraltypes.itp, and improper.itp into the ffbonded.itp doc in your CHARMM36 FF file. There are 4 main sections of the ffbonded.itp file; add bondtypes.itp to the [bondtypes] section, angletypes.itp to the [angletypes] section, dihedraltypes.itp to the first [dihedraltypes] section, and improper.itp to the second [dihedraltypes] section (it should mention impropers).
Note: In bondtypes.itp/angletypes.itp/etc., you only need to add the lines with information about peptoids. I have included the header for each section and followed each section with "; Original CHARMM..." so you know where to place the information we've gathered.

Lastly, add 'TC   12.01100 ; peptoid carbonyl carbon' to atomtypes.atp in your FF directory.

Once these updates have been made, you should be able to use the peptoid residues found in residue_pdb to build your own peptoid amino acid mimics and produce topologies using the updated CHARMM36 FF.
