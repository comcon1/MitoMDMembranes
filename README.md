# About

The repository contains data for modelling of slightly oxidized model mitochondrial membrane. Here we put models of monohydroperoxides of tetralinoleoyl cardiolipin as well as models of two BODIPY-based dyes sensitive to lipid peroxidation.

The deposited data supplements the original research:

>K. G. Lyamzaev, N. V. Sumbatyan, A. M. Nesterenko, E. G. Kholina, V. Natalia, S. Heinz-Jürgen, A. Y. Mulkidjanian, and B. V. Chernyak, **Mitoclox: A novel mitochondria-targeted fluorescent probe for tracing lipid peroxidation**, *Oxidative Medicine and Cellular Longevity*, vol. 2019, p. 9710208, 2019. [[doi](http://dx.doi.org/10.1155/2019/9710208)]

# Example system

This is an example system containing PLPC, PSLPE, POPC, oxidized TCLC and MitoCLOx dye.

![](https://github.com/comcon1/MitoMDMembranes/raw/master/model.bilayers/oxMembrane_hyp13_MitoCLOx/oxMembrane_hyp13_MitoCLOx_280ns.png)

# Repository structure
Into the `ff.construction` folder we put the protocols for parametrization for BODIPY-based dye and hydroperoxide of lynoleic acid. The `forcefield` folder contains ready GROMACS topologies and all files needed for compiling GROMACS binary file with the sample system. In the `model.bilayers` there are many subfolders with different types of mitochondrial model membranes, which differs with the type of lynoleic acid hydroperoxides in the TLCL molecules.


# Authors

**Alexey Nesterenko** 
* ResearcherID: [C-8646-2012](https://publons.com/researcher/C-8646-2012/)
* Scholar Profile: [A.M. Nesterenko](http://bit.ly/c1_scholar)
* RG profile: [Alexey Nesterenko](https://www.researchgate.net/profile/Alexey_Nesterenko)

**Ekaterina Kholina** 
* ResearcherID: [F-8720-2018](https://publons.com/researcher/F-8720-2018/)
