# sfxPhasing

MR_phasing stands for Molecular Replacement Automated Phasing Pipeline. It includes the Phenix.MRage and automatically for the grid search for user's defined parameter.

SAD_phasing stands for the Single-wavelength anomalous dispersion Automated Phasing Pipeline with heavy atoms of Selenium (Se) or Sulfur (S). It includes 2 separated modules: SHELXC/D, Crank2.
Later the Autobuild refinement will be considered to further refine the protein model.

Examples contain the SAD example of Streptavidin and the MR example of 4N5R which is a lysozyme example.

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
