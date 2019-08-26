from __future__ import print_function
import sys
import subprocess
import os
import argparse
import re
import json
import ast

if os.path.isfile("crank2.inp"):
    os.remove("crank2.inp")

parser= argparse.ArgumentParser()


parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-pdb","--pdb-file",help="input the pdb file",type = str)
parser.add_argument("-seq","--sequence-file",help="input the sequence file",type = str)


args = parser.parse_args()

if args.reflection_mtz:
    reflection_file = args.reflection_mtz
else:
    print('The reflection file is missing')

if args.pdb_file:
     pdb= args.pdb_file
else:
    print('The pdb file is missing')

if args.sequence_file:
    seq_file = args.sequence_file
else:
    print('The sequence file is missing')


os.system("cmtzsplit -mtzin "+reflection_file+" -mtzout reflection_out.mtz -colin 'F(+),SIGF(+),F(-),SIGF(-)' -colout Fplus,SIGFplus,Fminus,SIGFminus")

num_atoms = str(4)
with open('Guessed_atom_number.txt','r') as f:
    for line in f:
        if 'exp_num_atoms' in line:
	   num_atoms = line.split(" = ")[1].replace('\n','')
	   print(type(num_atoms))

line1 = 'fsigf plus dname=peak f=Fplus sigf=SIGFplus "file=reflection_out.mtz"'
line2 = 'fsigf minus dname=peak f=Fminus sigf=SIGFminus'
line3 = 'model substr atomtype=SE "file='+pdb+'" d_name=peak exp_num_atoms='+num_atoms+' monomers_asym=4 solvent_content=0.474'
line4 = 'sequence "file='+seq_file+'"'
line5 = 'handdet dmfull dm solomon phcomb multicomb'
line6 = 'dmfull threshold_stop::0.58 dm parrot phcomb refmac'

crank_inp_file = [line1, line2, line3, line4, line5, line6]

for i in crank_inp_file:
    print(i,file=open("crank2.inp", "a"))

command = 'python /reg/common/package/ccp4/ccp4-7.0/share/ccp4i2/pipelines/crank2/crank2/crank2.py --keyin crank2.inp --hklout result.mtz --xyzout result.pdb'

os.system(command)
