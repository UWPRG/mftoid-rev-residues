# mftoid-rev-residues

## Pfaendtner + Ferguson Research Team Collaboration

## MoSiC-CGenFF-NTOID: A Peptoid Simulation Package based on Weiser & Santiso's CGenFF-NTOID. Force Field Parameters, Structure Generator, and Simulation Scripts.

## REPOSITORY CONTENTS
	
1. merged.rtp: file containing atoms, bonds, atomtypes, and charges for each peptoid amino acid mimic
2. atomtypes.atp: file containing the new atomtypes needed for the augmented force field
3. ffnonbonded.itp: file containing the carbonyl carbon peptoid parameter from Mirijanian's MFTOID
4. ffbonded: directory containing 4 files (bondtypes.itp contains updated bond types, angletypes.itp angle types, dihedraltypes.itp proper dihedral types, and improper.itp the improper dihedral types)
5. residue_pdb: a directory containing PDBs of each of the peptoid versions of the amino acids, orderedand labeled according to updated forcefield residues found in merged.rtp
6. structure_maker: a directory used to generate structures for your simulated peptoids
7. simulation_template: a directory of files used to run MD simulations.
8. simulations: a directory in which your simulations will be stored if you choose to run our simulation scripts.
9. run_sim.sh: a script you can use to generate peptoid structures from sequence and run MD simulations in water.
10. environment.yml: an Anaconda environment you should use to run our scripts.
11. vacuum_sim.sh: a script to run simulations in vacuum
12. pyscf_sp: a folder containing our PySCF scripts

## INSTRUCTIONS

### Packages and Dependencies

The code in this package and in our scripts requires a set of dependencies. You can use our conda environment to have all of them properly installed. Once your github is cloned, this is done as follows: 
```
conda env create -f environment.yml
conda activate peptoid_env
```

### Setting up the augmented CHARMM36 Force Field

Download the CHARMM36 forcefield from the [MacKerell Lab website](https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs:~:text=charmm36%2Dfeb2021.ff.tgz). We highly recommend you use the February 2021 version of CHARMM36 to be compatible with our parameters and scripts into this Github. After downloading your CHARMM36 FF, the new parameters manually (as explained below) or automatically in the following way. First, ensure your GMXLIB environment variable is the directory that contains your CHARMM36-Feb2021 force field. Additionally, ensure you have Python3 installed. Second, ensure you have write permissions to the CHARMM36-Feb2021 force field files. Finally, enter the MFToid-Rev-Residues directory containing the contents of this github, and enter the following command: 

```
python attach_peptoid_ff.py
```

Otherwise, you may update your force field files manually:

1. Add the peptoid parameter TC found here in ffnonbonded.itp into the ffnonbonded.itp doc in your CHARMM36 FF file.

2. Add the information contained here in the ffbonded directory to the ffbonded.itp doc in your CHARMM36 FF file. There are 4 main sections of the ffbonded.itp file; add bondtypes.itp to the [bondtypes] section, angletypes.itp to the [angletypes] section, dihedraltypes.itp to the first [dihedraltypes] section, and improper.itp to the second [dihedraltypes] section (it should mention impropers).

Note: In bondtypes.itp/angletypes.itp/etc., you only need to add the lines with information about peptoids. I have included the header for each section and followed each section with "; Original CHARMM..." so you know where to place the information we've gathered.

3. Lastly, add 
```
TC   12.01100 ; peptoid carbonyl carbon
```
to atomtypes.atp in your FF directory.

Once these updates have been made, you should be able to use the peptoid residues found in residue_pdb to build your own peptoids and produce topologies using the updated CHARMM36 FF.

### Using the Structure Generator

Using the command prompt, enter the structure_maker directory. Type the following command, and a structural file will be created for your use (the name of which is customizable): 

```
python make_structure.py --seq [SEQUENCE] --mini [MINIMUM CODE] --file [FILENAME]
```
OR 
```
python make_structure.py --seq [SEQUENCE] --phi [PHI] --psi [PSI] --omega [OMEGA] --file [FILENAME]
```
Notes:
1. SEQUENCE: a string using the one-letter codes of any amino acid sequence. Additional residues exist, and their codes follow:  
    	+: positively charge histidine  
        @: histidine with a proton on the epsilon-nitrogen, rather than the delta-nitrogen  
   	X: aminoethyl (Nae) side chain, positively charged  
   	Z: S-1-phenylethyl (Nspe) side chain  
   	1: tert-butyl side chain (Ntbu)
   	2: phenyl side chain (Nph -- **DO NOT USE: UNSUITABLE**)  
   	3: p-bromophenyl-methyl side chain (Nbrpm)  
   	4: p-bromophenyl-ethyl side chain (Nbrpe)  
  	5: p-chlorophenyl-methyl side chain (Nclpm)  
   	6: p-chlorophenyl-ethyl side chain (Nclpe)  
   	7: p-fluorophenyl-methyl side chain (Nfpm)  
   	8: p-fluorophenyl-ethyl side chain (Nfpe)  
   	9: p-iodophenyl-methyl side chain (Nipm)  
   	0: p-iodophenyl-ethyl side chain (Nipe)
   	J: ethyl-methoxyethyl side chain (Neme)
   	U: decyl side chain (Ndec)
   	O: 1-hydroxyethyl sidechain (Noe)
3. MINIMUM CODE: one of the following twelve codes for the twelve possible minima:  
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
	HELIX: alpha-helix

	Please note not every amino acid sequence has a valid structure in any given minimum (proline, for example, may cause serious problems): be sure to use a valid minimum for your sequence.  
4. FILENAME: any filename is valid; however, if the filename does not end with ".pdb", that extension will be added automatically.

You may then use [FILENAME].pdb in your future simulations.

### Running Automated Simulation Scripts
We currently have simulation scripts set up to run 500 ns parallel bias metadynamics simulations in SPC/E water, biasing on all the backbone dihedral angles, radius of gyration, and alphabetas (distances) from the 12 common dihedral angle free energy minima. **These scripts generate the structure for you, so you do not have to run the structure maker code before running the simulation.** To run a simulation, use the following bash command. 
```
bash run_sim.sh [SEQUENCE] [MINIMUM CODE]
```
For practicality and speed concerns, we recommend you create an additional script to wrap this code so you may run it on a GPU or multiple parallel CPUs.   
To adjust conditions of your simulations, view the files in the ```simulation_template``` directory. Specifically, ```simulation_template/mini_sim.sh``` will be useful to modify for the number of molecules or the dimension of your box. The specific MD parameter files (.mdp) for your simulations are located in `2_em`, `3_md`, and `4_pb` in the `simulation_template` directory.

### Testing the High Throughput Simulation Package

To test the functionality of this simulation package, you may run the following bash command: ```bash test_example.sh```
This should run a simulation of disarcosine in water. After energy minimization and NVT and NPT equilibration, it will run for 5 picoseconds. It will create a folder called ```test_sim``` and the filenames inside that folder should be the same as the comparison folder ```test_sim_comp```. If you have any issues, please open an issue in this repository.
