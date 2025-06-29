## MoSiC-CGenFF-NTOID: A Peptoid Simulation Package based on Weiser & Santiso's CGenFF-NTOID. Force Field Parameters, Structure Generator, and Simulation Scripts. A Pfaendtner + Ferguson Research Team Collaboration   

View our preprint before using the force field: https://arxiv.org/abs/2409.06103.   

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
Additionally please use **GROMACS 2021.4** for the structure generation.

### Setting up the augmented CHARMM36 Force Field

Download the CHARMM36 forcefield from the [MacKerell Lab website](https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs:~:text=charmm36%2Dfeb2021.ff.tgz). We highly recommend you use the February 2021 version of CHARMM36 to be compatible with our parameters and scripts into this Github. After downloading your CHARMM36 FF, the new parameters manually (as explained below) or automatically in the following way. First, ensure your GMXLIB environment variable is the directory that contains your CHARMM36-Feb2021 force field. Additionally, ensure you have Python3 installed. Second, ensure you have write permissions to the CHARMM36-Feb2021 force field files. Finally, enter the MFToid-Rev-Residues directory containing the contents of this github, and enter the following command: 

```
python attach_peptoid_ff.py
```
For the sake of the compatibility of GROMACS and CHARMM36 (see: https://gromacs.bioexcel.eu/t/capping-residue-ace-not-recognised-from-the-rtp-file/1694), replace the `[ ACE ]` residue (**NOT ACEp**) in your `merged.rtp` with the following text:
```
[ ACE ] ; N-terminal acetyl patch
 [ atoms ]
           CAY     CT3     -0.270  1
           HY1     HA3     0.090   1
           HY2     HA3     0.090   1
           HY3     HA3     0.090   1
           C       C       0.510   2
           O       O       -0.510  2
 [ bonds ]
           C       CAY
           C       +N
           CAY     HY1
           CAY     HY2
           CAY     HY3
           O       C
 [ impropers ]
            C      CAY     +N      O
```

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

### Adding new residues

You may want expand this package to accommodate simulations of peptoids of interest to you, which may include residues not yet incorporated into this package. These resiudes may be added by first adding the necessary force field parameters in the same way as one would do with general CHARMM, and then (optionally) adding files and codes enabling the compatibility with the structure generator. First, for the force field paramters:   
1. Create a residue topology (.rtp) file for your residue, and append it to merged.rtp
2. Add any atom types not familiar to CHARMM36-Feb2021 into atomtypes.atp   
3. Use the CGenFF parameter generator to determine which additional bonded force field parameters need to be added. Add the resulting bond, angle, and dihedral parameters, **with the exception of phi, psi, omega, and rho dihedral angles, as well as chi_1 dihedral angles with atom types identical to NPhe or Nspe**, to ffbonded.itp in the same locations as the other MoSiC-CGenFF-NTOID parameters. **Note that you must convert the force constants from CHARMM to GROMACS-compatible units**.
4. Add a PDB of your residue into the folder "residue_pdb" with the the filename [RESp].pdb, where [RESp] is replaced by the name of your residue.   

Then, optionally for the structure generator:
1. Copy the PDB file of your residues and remove the backbone atoms and alpha-hydrogens of the copied PDB. Save it in the structure_maker/.nobkb_residue_pdb folder as [RESp].pdb, where [RESp] is replaced by the name of your residue.
2. Create a new single-character code for your residue (for a single residue, we suggest 'B'. Beyond this, you may have to get creative with special characters or lowercase letters.) Add a line to the FILENAMES dictionary of make_structure.py  with the following information: '[character]: [RESp]', where [RESp] is replaced by the name of your residue.

You may now take full advantage of our high-throughput simulation package for your new residue.

### Testing the High Throughput Simulation Package

To test the functionality of this simulation package, you may run the following bash command: ```bash test_example.sh```
This should run a simulation of disarcosine in water. After energy minimization and NVT and NPT equilibration, it will run for 5 picoseconds. It will create a folder called ```test_sim``` and the filenames inside that folder should be the same as the comparison folder ```test_sim_comp```. If you have any issues, please open an issue in this repository.

### Citation
```
@misc{berlaga2024,
      title={A Modular and Extensible CHARMM-Compatible Model for All-Atom Simulation of Polypeptoids}, 
      author={Alex Berlaga and Kaylyn Torkelson and Aniruddha Seal and Jim Pfaendtner and Andrew L. Ferguson},
      year={2024},
      eprint={2409.06103},
      archivePrefix={arXiv},
      primaryClass={cond-mat.soft},
      url={https://arxiv.org/abs/2409.06103}, 
}
```
### Upcoming updates
Incorporating the ACE changes within `attach_peptoid_ff.py`, adding hydrogen database (.hdb) entries for `pdb2gmx` using hydrogen-free input files.
