# sfxPhasing

MR_phasing refers to Molecular Replacement Automated Phasing Pipeline. It runs Phenix.MRage with automatic grid search and/or with user defined parameters.

SAD_phasing refers to Single-wavelength Anomalous Dispersion Automated Phasing Pipeline with heavy atoms of Selenium (Se) or Sulfur (S). It uses 2 programs: SHELXC/D, Crank2.
Later Autobuild refinement may further refine the protein model.

Examples contain a Se-SAD example of Selenobiotinyl-Streptavidin and an MR example of 4N5R lysozyme.

## Required External Packages:
The following packages are required to run sfxPhasing scripts.

phenix version: 1.16
Freely available after requesting a download password: https://www.phenix-online.org/download/

ccp4 version: 7.0.077 (It includes SHELXC version:2016/1, SHELXD version:2013/2, and Crank2 version:2.0.229)
Freely available: http://www.ccp4.ac.uk/download

pymol version: 2.3.4
Freely available: https://pymol.org/2/

## Python version:
python: 2.7

## Datasets:

Download example data sfxPhasing.zip from here (by invitation only): 
https://stanford.box.com/s/8s55n5smgxb7gt00ihzoj241qgk1bocl

Uncompress the data folder:  
$unzip sfxPhasing.zip
