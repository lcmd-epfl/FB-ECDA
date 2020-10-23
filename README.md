# FB-ECDA [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4121975.svg)](https://doi.org/10.5281/zenodo.4121975)

## Introduction
**FB-ECDA** is a decomposition analysis tool for electronic coupling term in charge transfer process of organic semiconductors. Please cite
our work if you find this tool is useful (the work is to be submitted soon).The current version supports computations of 
fragment-fragment electronic coupling between two molecules of the same kind. It also supports computations of electronic
coupling of arbitrary MO by specifying its index. If one is interested in hole (electron) transfer, please specify the index
of the HOMO (LUMO) in Gaussian log file.

## Usage
To compute the fragment-fragment electronic coupling, one needs to first perform single-point computations for the dimer
and two monomers in the dimer (the coordinate should be the same, do not rotate the molecules). The Gaussian input setting
can be found in **example.tar.gz**.

The usage can be found by executing the python file FB-ECDA.py with -h:    
`python FB-ECDA.py -h`  

In the example file, one can find the Gaussian input (**.com**) that is recommended for performing FB-ECDA.
In addition, you can find the input necessary for performing FB-ECDA computation:

|Input files|Explanation|
| :---: | --- |
|13.log|Gaussiagn log file. The overlap matrix and Fock matrix of molecule 1.|
|169.log|Gaussiagn log file. The overlap matrix and Fock matrix of molecule 2.|
|13-169.log|Gaussiagn log file. The overlap matrix and Fock matrix of molecule 1+2.|
|13.pun|Gaussian punch file. Coefficient and eigenvalue of each MO of molecule 1.|
|169.pun|Gaussian punch file. Coefficient and eigenvalue of each MO of molecule 2.|
|13-169.pun|Gaussian punach file. Coefficient and eigenvalue of each MO of molecule 1+2.|
|AO.txt|User specified Input. Number of atomic orbitals of each element for the chosen basis set.|
|frag.txt|User specified Input. Definition of fragments by specifying the index of atoms as specified in the Gaussian monomer input files.|

With these input files, we can perform FB-ECDA following a straightforward manner:          
`python FB-ECDA.py 13.log 169.log 13-169.log 13.pun 169.pun 13-169.pun 152 AO.txt frag.txt`     
where 152 is the index of HOMO for the molecule in the input file    

The output shows each fragment-fragment electronic coupling term (V_frag1-frag2) and the total electronic coupling (V_tot) in eV:

`V_core-core: 0.00043440836750638194`  
`V_core-arm: -0.0001615526719718881`  
`V_arm-core: 0.0006802747262879584`  
`V_arm-arm: -0.004705632444962128`  
`V_tot: -0.003752502023139676`  

## References

> To be submitted soon

## Support and development
For bug reports/suggestions/complaints please raise an issue on [GitHub].

Or contact us directly: [imagine7801@gmail.com]

[GitHub]:<https://github.com/pumachu/FB-ECDA>
[imagine7801@gmail.com]:<mailto:imagine7801@gmail.com>
