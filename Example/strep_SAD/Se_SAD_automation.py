from __future__ import print_function
import sys
import subprocess
import os
import argparse
import re
import json
import ast
import numpy as np
import subprocess
import time
import shutil

current_path = os.getcwd()

if os.path.isfile("all_jobs.txt"):
    os.remove("all_jobs.txt")
    
parser= argparse.ArgumentParser()

parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-seq","--sequence-file",help = "input the sequence file", type = str)
parser.add_argument("-FIND", "--number-of-atoms", default = 4, help='enter the number of atoms', type = str)
parser.add_argument("-resl","--resolution", default = '3.0', help = "input resolution", type = str)

args = parser.parse_args()

if args.reflection_mtz:
   reflectionFile = args.reflection_mtz
   job_name = args.reflection_mtz.replace('.mtz','')
else:
   print('Please input the mtz file for SHELXC/D')
   sys.exit()

if args.sequence_file:
   sequenceFile = args.sequence_file
else:
   print('Please input sequence for crank2')
   sys.exit()

if args.resolution:
   resolution = args.resolution

if args.number_of_atoms:
   atom_find = args.number_of_atoms

shelx_CD = 'python SHELX_script.py  -rfl '+reflectionFile+' -MIND1 -3.5 -MIND2 2.2 -resl '+str(resolution)+' -TEST 0 99 -ESEL 1.5 -SFAC SE -NTRY 1000 -FIND '+str(atom_find)

os.system(shelx_CD)

# To check whether SHELXC/D
if os.path.isfile(job_name+'_fa_cleaned.pdb'):
    new_pdb = job_name+'_fa_cleaned.pdb'
else:
    print('cleaned pdb has not been properly generated')
    os.chdir('../')
    shutil.rmtree(current_path)
    sys.exit()

# run crank2
crank2_cl = 'python crank2_script.py -rfl '+reflectionFile+' -pdb '+new_pdb+' -seq '+sequenceFile
os.system(crank2_cl)
# check whether crank2 works
os.chdir('crank2')
content = open('crank.loggraph','r').read()
for line  in content.split('\n'):
    if 'The final FOM is' in line:
        FOM = float(line.split('is ')[1])
        if FOM < 0.4:
           os.chdir(current_path)
           os.chdir('..')
	   shutil.rmtree(current_path)
	   sys.exit()
        else:
	   os.chdir(current_path)
           print(current_path+'/'+str(FOM),file = open("FOM.txt","a"))
           if os.path.isfile('result.mtz'):
	      new_mtz = 'result.mtz'
	   else:
	      print("experimental reflection file has not been genereted")
	      sys.exit()

auto_build_cl = 'python autobuild.py -rfl '+new_mtz+' -seq '+sequenceFile+' -rfff 0.05 -nproc 12'

os.system(auto_build_cl)

os.chdir('AutoBuild_run_2_/')

list=[]
with open('AutoBuild_summary.dat','r') as f:
    for line in f:
        list.append(line.rstrip('\n'))

os.chdir('..')
for i in list:
    if 'Best solution on cycle' in i:
        print(os.getcwd()+'/'+i, file = open("best_solution_R_Rfree.txt",'a'))

#run autobuild
# if os.path.isfile('result.mtz'):
#     new_mtz = 'result.mtz'
# else:
#     print('experimental reflection file has not been properly generated') 
#     sys.exit()

# auto_build_cl = 'python autobuild.py -rfl '+new_mtz+' -seq '+sequenceFile+' -rfff 0.05 -nproc 4'

# os.system(auto_build_cl)



