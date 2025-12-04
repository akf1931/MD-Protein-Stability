# MD-Protein-Stability

## Summary

Molecular dynamics (MD) simulations of Green Fluorescent Protein (GFP) at baseline condition and stress condition with analysis.

## Procedure

### Step 1: A simulation-ready .pdb file of GFP

Download the 1EMA.pdb file.
Open the file with Notepad++.
In Notepad++, use ctrl+F to open the 'Find' window and search for 'MISSING RESIDUES' to locate the following lines:

REMARK 465   M RES C SSSEQI                                                     
REMARK 465     MSE A     1                                                      
REMARK 465     THR A   230                                                      
REMARK 465     HIS A   231                                                      
REMARK 465     GLY A   232                                                      
REMARK 465     MSE A   233                                                      
REMARK 465     ASP A   234                                                      
REMARK 465     GLU A   235                                                      
REMARK 465     LEU A   236                                                      
REMARK 465     TYR A   237                                                      
REMARK 465     LYS A   238                                                      

These REMARKs indicate the protein is missing its first residue and the entire 230 - 238 loop.

Beneath the MISSING RESIDUES section is the MISSING ATOM section:

REMARK 470   M RES CSSEQI  ATOMS                                                
REMARK 470     GLU A   6    CD   OE1  OE2                                       
REMARK 470     LYS A  26    NZ                                                  
REMARK 470     LYS A  52    CD   CE   NZ                                        
REMARK 470     LYS A 101    NZ                                                  
REMARK 470     LYS A 107    CE   NZ                                             
REMARK 470     ARG A 122    NE   CZ   NH1  NH2                                  
REMARK 470     GLU A 124    CD   OE1  OE2                                       
REMARK 470     LYS A 131    CE   NZ                                             
REMARK 470     GLU A 132    CG   CD   OE1  OE2                                  
REMARK 470     LYS A 156    CG   CD   CE   NZ                                   
REMARK 470     GLN A 157    CB   CG   CD   OE1  NE2                             
REMARK 470     LYS A 158    CB   CG   CD   CE   NZ                              
REMARK 470     LYS A 162    CE   NZ                                             
REMARK 470     ASN A 212    OD1  ND2                                            
REMARK 470     LYS A 214    CE   NZ 

These REMARKs indicate the protein is missing sidechain atoms for 15 residues. We know they are sidechain atoms for two reasons. The first is that none of the ATOMS are C, O, CA, or N, which are the backbone atom names used in PDB files. Second, the authors of the crystal structure tell us so in the header for the file:

REMARK   3  OTHER REFINEMENT REMARKS:                                           
REMARK   3  THE FINAL (FO-FC) DENSITY SHOWS LARGE DIFFERENCE FEATURES           
REMARK   3  LOCATED AROUND THE MAIN CHAIN PART OF THE FLUOROPHORE               
REMARK   3  THAT ORIGINATES FROM RESIDUES THR 65 AND GLY 67.  THE               
REMARK   3  DIFFERENCE DENSITIES ALSO AFFECT THE ATOMIC-POSITION                
REMARK   3  REFINEMENT OF VAL 68.  THE DIFFERENCE FEATURES MIGHT BE             
REMARK   3  EXPLAINED AS A SUBSET (<30%) OF MOLECULES THAT HAS FAILED           
REMARK   3  TO UNDERGO THE COMPLETE FORMATION OF THE FLUOROPHORE.  THE          
REMARK   3  LOOP RESIDUES 157 AND 158 ARE DISORDERED.  A NUMBER OF              
REMARK   3  SURFACE RESIDUES HAVE TRUNCATED SIDE CHAINS DUE TO WEAK             
REMARK   3  OR NO DENSITY.                         

To prepare the PDB file for simulation, coordinates for the missing residues and atoms must be supplied. There are many strategies and tools helpful for accomplishing this task. For this experiment, access the Google Colab AlphaFold2.ipynb notebook (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) to model the entire *A. victoria* GFP sequence, obtained from the .fasta file provided at the 1EMA PDB repository.

Open the 1EMA.pdb and the AlphaFold2 rank 1 structure files in PyMOL. 
Align the rank 1 output structure with the 1EMA structure using PyMOL's plugin alignment wizard, using a one-to-one alignment with default options. Copy the AlphaFold2 structure to a new object (obj01).
Select the three residues the unusual fluorophore CRO residue replaces in obj01: 65, 66, and 67. Use an action on the selection to delete the atoms. Select the CRO residue in the 1EMA structure and copy the selection to obj01. This will append the CRO residue to the end of the AlphaFold2 structure sequence. Select it there. In the PyMOL command interface, run the commands:

	alter sele, segi=""
	alter sele, chain="A"

This will remove the segment identifier PyMOl automatically assigned to the copied CRO residue and change its chain identifier to A, to match the rest of the residues in obj01.

Copy obj1 to a new object (obj02). This will insert the CRO residue into its correct place in chain A via its residue ID number. 

To obtain a forcefield capable of simulating the unusual CRO residue that represents the fluorophore, download the charmm36-jul2022 forcefield for GROMACS: https://mackerell.umaryland.edu/charmm_ff.shtml. Place it in the gromacs/share/top folder alongside the other default .ff forcefield folders.

Now the atom names of the CRO residue must match the [ CRO ] entry in the charm37-jul2022 aminoacids.rtp file. This is a bit like playing a puzzle game, as the only structural information provided by the [ CRO ] entry are the atom bonds. One atom at a time, ignoring the absent hydrogen atoms for now, use the PyMOL command with only the target atom selected:

	alter sele, name="XXX"

Replacing XXX with the appropriate atom name from the [ CRO ] entry.

Once all atom names are assigned and logically sound with the bond information from the [ CRO ] entry, use the PyMOL builder tool to add hydrogens. Unfortunately, CRO does not have an entry in any .hdb file, so the hydrogens must be added manually and named appropriately to what the entry in the .rtp file expects. A quick way to do this is to use the builder tool's 'Add H' button, which will add hydrogens to the entire obj02. Then run the command

	select to_delete, element H and not resn CRO

Then use an action on the selection to delete the selected atoms. These hydrogens should be removed because pdb2gmx will supply them later, guaranteeing the use of the expected hydrogen names for the forcefield.

Four more hydrogens must be deleted. Two are the extra hydrogen on CRO's 'N' atom, leaving it with 1, and the hydrogen on CRO's 'C' atom, leaving it with no hydrogens. These hydrogens do not exist in the [ CRO ] .rtp file entry because the bonds they take up are meant to be occupied with continuing the backbone of the protein. The third hydrogen to delete is on CRO's 'O2' atom, as this is a carbonyl oxygen participating in the cyclic amide. The fourth hydrogen is on the 'CB2' atom, as it is conjugated alkene and requires only one hydrogen.

Now use the alter sele command to name the hydrogens in accordance to their bond information in the [ CRO ] entry of the .rtp file.

Export the grafted object, obj02, as the GFP_A_victoria.pdb file. Delete the TER record PyMOL's export function automatically adds around the CRO residue lines, as this will confuse pdb2gmx into thinking there are two chains if left in.

The GFP_A_victoria.pdb structure is simulation-ready, with no crystal waters, with all crystallographic MSE residues replaced with the more biologically representative MET residues, and with the fluorescent CRO residue affixed.

To prepare the file structure for conducting molecular dynamics simulations for behavioral comparison of GFP in physiological and stress conditions, create a data/300 K 150 mM and a data/340 K 1000 mM folder. Within each of those folders, make r-1, r-2, and r-3 folders to represent three replicates of each condition.

### Step 2: A 200 ns molecular dynamics trajectory under physiological conditions

Add the line 'CRO Protein' to the top of the residuetypes.dat file in the /gromacs/share/top folder so pdb2gmx will know to treat it as part of the protein structure. In \data\300 K 150 mM\r-1, run the command to prepare the topol.top topology file, posre.itp position restraint file, and the .gro structure file from a copy of the prepared GFP_A_victoria.pdb file.

	gmx pdb2gmx -f GFP_A_victoria.pdb -o processed.gro -water tip3p -ter

Select the charmm36-jul2022 forcefield and designate both termini to be charged. The command output tells us the total atomic mass of the protein is 26,902.531 a.m.u., and the total charge of the protein is -6 e.

We now define the periodic boundary conditions for the system. These define the size of the tiny universe the simulation models. The protein will drift freely in water (still to be added), and when any of its atoms cross a periodic boundary, they will reappear on the other size of the universe. In this way, we do not impose artificial walls the protein can encounter while still putting a limit on the number of atoms we simulate.

	gmx editconf -f processed.gro -o newbox.gro -c -bt dodecahedron -d 1.0

The above command centers the protein in the simulation box (-c), defines the boundaries as forming a dodecahedron shaped universe- which is an efficient shape for globular proteins (-bt dodecahedron)- and ensures at least 10 angstrom of space are between the boundaries and the outermost atoms of the protein (-d 1.0), so it cannot encounter itself as it drifts across the boundaries.

Next, add solvent to the system- pure water using the tip3p water molecule definition, which is a good compromise between realism and processing efficiency.

	gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

The next command, genion, requires a .tpr file as input. To create one, supply the minim.mdp file with a maxwarn flag since it will not like the unbalanced charge of the system. As noted at https://manual.gromacs.org/2024.0/user-guide/force-fields.html, make sure to use the correct mdp file options when using the CHARMM forcefield.

	gmx grompp -f minim.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1

Neutralize the system and set 150 mM ionic strength, a physiologically comfortable salt content for GFP.

	gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

Select group 13, SOL, to replace some water molecules with ions.

Reuse the minim.mdp file in creating the energy minimization .tpr file.

	gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
	gmx mdrun -v -deffnm em

Apply the physiological nvt.mdp file in performing temperature equilibration: constant number of atoms (N), system volume (V), and temperature (T) at 300 K- a physiologically comfortable temperature for GFP.

	gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
	gmx mdrun -v -deffnm nvt

Apply the physiological npt.mdp file in performing pressure equilibration: constant number of atoms (N), pressure (P), and temperature (T), allowing volume to scale evenly to achieve 1 bar of pressure.

	gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
	gmx mdrun -v -deffnm npt

Apply the physiological md.mdp file to perform a 200 ns production run of molecular dynamics.

	gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
	nohup gmx mdrun -v -deffnm md &

### Step 3: Perform two more replicates of GFP under physiological conditions

Beginning from the energy minimized structure but in new folders, re-conduct physiological NVT, NPT, and MD.

	gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
	gmx mdrun -v -deffnm nvt
	gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
	gmx mdrun -v -deffnm npt
	gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
	nohup gmx mdrun -v -deffnm md &

### Step 4: A system prepared for molecular dynamics at stress conditions

Repeat steps 2 and 3 in a fresh folder using the stress .mdp files to perform the same 200 ns simulation now under stress conditions of 340 K and 1 M NaCl.

Starting with the solvated system, run the genion command to introduct 1 M salt.

	gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 1
	
Continue through EM and three repeats of NVT, NPT, and MD. The stressful 340 K temperature setting is defined in the .mdp files.

### Step 5: Clean the trajectories

Per trajectory, run the command to only retain every fifth frame while also centering the protein in the simulation box and discarding the waters and ions.

	gmx trjconv -f md.xtc -s md.tpr -o md_molcompact.xtc -ur compact -pbc mol -skip 5

Choose group 1, Protein, in the option menu presented to output only the protein structure.

To make sure the protein is not split across the periodic boundaries, we must use an example structure of it when it is whole. Use PyMOL to evalute identify a .gro file- which are output at each EM, NVT, and NPT step- that shows the GFP not split across any PBCs. The accompanying .tpr file to that .gro file can be used to make sure the GFP stays whole throughout the cleaned trajectory. First, we remove the waters and ions from the select .tpr file as well so it matches our md_molcompact.xtc trajectory.

	gmx convert-tpr -s npt.tpr -o npt_1.tpr

Choose group 1, Protein.

We now use the selected, edted .tpr file to run the -pbc nojump command of trjconv on the md_molcompact.xtc trajectory.

	gmx trjconv -s npt_1.tpr -f md_molcompact.xtc -pbc nojump -o md_molcompact-nojump.xtc

Chose group 0, System, as output- because now the system already is composed of only the protein.

Use the edited .tpr file again to perform a progressive least-squares structure fit frame by frame, so the protein remains centered in view when the trajectory is loaded into PyMOL. Note this command does not change the flexing of the protein, only its position in space. The output trajectory can therefore be used for statistical analysis on the protein's behavior over the trajectory and between different trajectories.

	gmx trjconv -s npt_1.tpr -f md_molcompact-nojump.xtc -fit progressive -o md_molcompact-nojump-progressive.xtc

Chose group 1, Protein, for least squares fit and group 0, System, for output.

Finally, we produce a .gro file to visualize and load the trajectory onto in PyMOL.

	gmx trjconv -s npt_1.tpr -f md_molcompact-nojump-progressive.xtc -dump 1000 -o md_1ns.gro

Chose group 0, System, for output.

### Step 6: Output measurements

From the /data folder, run Python script at scripts/extract_observables.py to record .csv files containing, per run, the measurements of radius of gyration (Rg), root mean square deviation (RMSD), helix secondary structure fraction, sheet secondary structure fraction, and other secondary structure fraction.

### Step 7: Analysis

Refer to the notebook/analysis.ipynb file for evaluating the raw measurements of each trajectory obtained above.
