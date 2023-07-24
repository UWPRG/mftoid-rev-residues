# mftoid-rev-residues

## Pfaendtner + Ferguson Research Team Collaboration

## GOAL: Build library of residues based on the Weiser+Santiso modified MFTOID force field peptoid backbone

### REPOSITORY CONTENTS
	
1. merged.rtp: file containing atoms, bonds, atomtypes, and charges for each peptoid amino acid mimic

2. ffbonded: directory containing 4 files (bondtypes.itp contains updated bond types, angletypes.itp angle types, dihedraltypes.itp proper dihedral types, and improper.itp the improper dihedral types)

3. ffnonbonded.itp: file containing the carbonyl carbon peptoid parameter from Mirijanian's MFTOID

4. residue_pdb: a directory containing PDBs of each of the peptoid versions of the amino acids, orderedand labeled according to updated forcefield residues found in merged.rtp

### INSTRUCTIONS

Download the CHARMM36 forcefield from the MacKerell Lab website; there have been some updates between February 2021 and July 2022 on forcefield file structure. We have been using the February 2021 version of the FF (and older) so we will most closely match that directory structure.

To update the downloaded version of your CHARMM FF for peptoid amino acid mimic capabilities, start by copying the information found in this Github repository in the merged.rtp file directly into the merged.rtp doc in your CHARMM36 FF file.

After downloading your CHARMM36 FF, the new parameters manually (as explained below) or automatically in the following way. First, ensure your GMXLIB environment variable is the directory that contains your CHARMM36 force field. Additionally, ensure you have Python3 installed. Then, enter the MFToid-Rev-Residues directory containing the contents of this github, and enter the following command: 

**python attach_peptoid_ff.py**

Otherwise, you may update your force field files manually:

Next, add the peptoid parameter TC found here in ffnonbonded.itp into the ffnonbonded.itp doc in your CHARMM36 FF file.

Next, add the information contained here in the ffbonded directory to the ffbonded.itp doc in your CHARMM36 FF file. There are 4 main sections of the ffbonded.itp file; add bondtypes.itp to the [bondtypes] section, angletypes.itp to the [angletypes] section, dihedraltypes.itp to the first [dihedraltypes] section, and improper.itp to the second [dihedraltypes] section (it should mention impropers).

Note: In bondtypes.itp/angletypes.itp/etc., you only need to add the lines with information about peptoids. I have included the header for each section and followed each section with "; Original CHARMM..." so you know where to place the information we've gathered.

Lastly, add 'TC   12.01100 ; peptoid carbonyl carbon' to atomtypes.atp in your FF directory.

Once these updates have been made, you should be able to use the peptoid residues found in residue_pdb to build your own peptoid amino acid mimics and produce topologies using the updated CHARMM36 FF.


### STRUCTURE MAKER USAGE INSTRUCTIONS

The structure maker creates a "pdb" file of a peptoid with a natural amino acid sequence in a desired minimum configuration. It makes use of Python 3.7, MDTraj 1.9.4, and Bio.PDB (biopython 1.79). The packages may be installed with Anaconda. 

Clone the repository locally or download the "structure_maker" directory. Then, using the command prompt, enter the structure_maker directory. Type the following command, and a structural file will be created for your use (the name of which is customizable): 

**python minima.py --seq [SEQUENCE] --mini [MINIMUM CODE] --file [FILENAME]**

Notes:
1. SEQUENCE: a string using the one-letter codes of any amino acid sequence. Additional residues exist for a positively charge histidine (code: B), histidine with a proton on the epsilon-nitrogen, rather than the delta-nitrogen (code: J), as well as a tert-butyl side chain (code: 1) and a phenyl side chain (code: 2).
2. MINIMUM CODE: one of the following twelve codes for the twelve possible minima:
  	A-T: alpha minus, trans
   	A-C: alpha minus, cis
   	A+T: alpha plus, trans
   	A+C: alpha plus, cis
    	AD-T: alpha-D minus, trans
   	AD-C: alpha-D minus, cis
	AD+T: alpha-D plus, trans
   	AD+C: alpha-D plus, cis
   	C7B-T: C7-beta minus, trans
   	C7B-C: C7-beta minus, cis
	C7B+T: C7-beta plus, trans
   	C7B+C: C7-beta plus, cis
   Please note not every amino acid sequence has a valid structure in any given minimum (proline, for example, may cause serious problems): be sure to use a valid minimum for your sequence.
3. FILENAME: any filename is valid; however, if the filename does not end with ".pdb", that extension will be added automatically.

You may then use [FILENAME].pdb in your future simulations.
