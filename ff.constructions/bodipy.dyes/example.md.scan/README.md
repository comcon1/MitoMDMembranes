# PES scan in MD

*ProfileMD* script builds energy profile as function of the dihedral angle
between the atoms of the molecule. It executes energy minimization and
constrained dynamics for the range of angle values and writes down resulting
average energies. 

## Authors

>initially written by  
**Tatiana V. Galochkina**  
email: <tat.galochkina@gmail.com>

>Université Sorbonne Paris Cité, Université Paris Diderot, Inserm, INTS  
Unité Biologie Intégrée du Globule Rouge UMR S1134, DSIMB, Laboratoire d’Excellence GR-Ex

## How to run?

To run this script you should have to fit the following requirements:

1. File `start.gro` should exist

    This is a regular structure file of your molecule. 

2. File `dih.ndx` should exist

    This is an index file with exactly one group containing exactly four atoms.
These four atoms define the angle of rotation. All other lines of the file are
ignored.

3. File `temp.top` should exist

    Your topology file should have `"#define mydih"` definition of improper dihedral 
at 10th line. In `[dihedrals]` section you should set to `"mydih"` your angle of
interest (the same as prescribed in your dih.ndx file). NOTE that all other
diherdrals around the same valence bond must be turned off. Please take a look
on lines 183-188 of the example file.


4. GROMACS version > 5.*

    Our script is designed for running all the commands through *gmx* application
wrapper.


5. You should have protocols of energy minimization and constrained dynamics.

    You can use *em.mdp* and *run.mdp* provided in our example. NOTE that you should not use PME to calculate electrostatic interactions for small system in a big box. It would lead the useless increase of the CPU time.


Once all the listed files are prepared you can run *ProfileMD*. To check if the
script gives a correct result, take a look on the output \_\*.pdb file. Please be
careful while preparing your topology and selecting terms from the energy
output. By default the script collects first 14 terms from EDR files. If you
include some special energy groups you may be interested in increasing this
number. In this case please address *GetEnergy* function.

All intermediate calculations are stored in the scratch folder. To repeat the
scan procedure you should remove scratch folder or change its name.

# How to try

Load force field from this repository, save it to `charmm36.ff` folder here and then run `ProfileMD.py` from the current directory. 
