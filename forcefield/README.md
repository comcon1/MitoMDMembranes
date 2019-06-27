We put the force field as well as the minimal example for building TPR-file.

* `toppar` -- included ITP files
     * TLXn -- monohydroperoxide of TLCL. @todo: link to folder with RTP
     * mitoperox -- topology of MitoPerOx (Prime et al., 2012: [doi](https://doi.org/10.1016/j.freeradbiomed.2012.05.033))
     * mitoperox2 -- topology of MitoCLOx (in press)

* `patch-for-charmm36` -- Patch for CHARMM36 folder
   * `ffbonded-c1.itp` -- additional bonded interactions 
   * `atomtypes-c1.itp` -- additional atomtypes
   * `ffnonbonded-c1.itp` -- additional non-bonded interactions
   * `forcefield.itp` -- rewrite of main ff file: *ffbonded-c1* and *ffnonbonded-c1* are included. 

* `itpFromCHARMMGUI` -- lipid topologies used by us and downloaded from CHARMM-GUI are deposited in this folder. Following lipid topologies were used:
   * TLCL2 -- tetralinoleoylcardiolipin (charge -2)
   * PLPC -- 1-palmitoyl-2-linoleoyl-phosphocholine
   * POPC -- 1-palmitoyl-2-oleoyl-phosphocholine
   * SLPE --1-stearyl-2-linoleoyl-phosphoethanolamine

* `membrane.gro` -- membrane example
* `topol.top` -- topology for example membrane

Please download [original force field at MacKerell lab site](http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs) and put patch files into the folder. We do not put the force field here because of license limitations. One can use any version of CHARMM forcefield derived from *CHARMM36-nov2016* if these files are added into. Also download the following lipid files from CHARMM-GUI into `toppar` folder.

Please run the following command to check the force field:

`gmx grompp -c membrane.gro -n index.ndx -f protocol.mdp -p topol.top -o t.tpr`

If the command finalizes without errors then FF is correct.
