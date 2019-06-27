We put the force field as well as the minimal example for building TPR-file.

* `toppar` -- included ITP files
   * Lipids from CHARMM-GUI:
     * TLCL2 -- tetralinoleoylcardiolipin (charge -2)
     * PLPC -- 1-palmitoyl-2-linoleoyl-phosphocholine
     * POPC -- 1-palmitoyl-2-oleoyl-phosphocholine
     * SLPE --1-stearyl-2-linoleoyl-phosphoethanolamine
   * Our molecules:
     * TLXn -- monohydroperoxide of TLCL. @todo: link to folder with RTP
     * mitoperox -- topology of MitoPerOx (Prime et al., 2012: [doi](https://doi.org/10.1016/j.freeradbiomed.2012.05.033))
     * mitoperox2 -- topology of MitoCLOx (in press)
* `charmm36-mar2019.ff` -- [Original force field at MacKerell lab site](http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs)
   * `ffbonded-c1.itp` -- additional bonded interactions 
   * `atomtypes-c1.itp` -- additional atomtypes
   * `ffnonbonded-c1.itp` -- additional non-bonded interactions
   * `forcefield.itp` -- *ffbonded-c1* and *ffnonbonded-c1* are included. 

One can use any version of CHARMM forcefield derived from *CHARMM36-nov2016* if these files are added into.

* `membrane.gro` -- membrane example
* `topol.top` -- topology for example membrane

Please run the following command to check the force field:

`gmx grompp -c membrane.gro -n index.ndx -f protocol.mdp -p topol.top -o t.tpr`

If the command finalizes without errors then FF is correct.
