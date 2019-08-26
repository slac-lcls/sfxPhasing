import sys
import shutil
import subprocess
import os
import argparse
import re
import numpy as np

parser= argparse.ArgumentParser()

refl_file_input = ''
pdb_file_input = ''
seq_file_input = ''
data_labels = ''

parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-resl","--resolution", default = "2.5", help="define the resolution", type = str)
parser.add_argument("-pdb","--pdb-file",help="input the pdb file",type = str)
parser.add_argument("-seq","--sequence-file",help="input the sequence file",type = str)

args = parser.parse_args()

if args.reflection_mtz:
    refl_file_input = args.reflection_mtz
else:
    print('The reflection file is missing')

if args.pdb_file:
     pdb_file_input= args.pdb_file
else:
    print('The pdb file is missing')

if args.sequence_file:
    seq_file_input = args.sequence_file
else:
    print('The sequence file is missing')

process = subprocess.Popen('phenix.xtriage '+refl_file_input+' '+pdb_file_input+' '+seq_file_input, 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()

split_out=out.splitlines()
my_list=[]
for i in range(len(split_out)):
    my_list.append(split_out[i].decode("utf-8"))
    #print(split_out[i].decode("utf-8"))

for line in my_list:
    if 'Best guess' in line:
        string = line
        matthew_coefficient = int(re.findall('\d+', string )[0])

    if "Data labels" in line:
        data_labels = line.split(':')[1].replace(' ','')

#Change the range in the list below
rmsd = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
#rmsd = [1.0,1.1]

if matthew_coefficient == 1:
    mat_coef_grid = range(1,3)
else:
    mat_coef_grid = range(matthew_coefficient-1,matthew_coefficient+2)

for directory in os.listdir(os.getcwd()):
    if 'Matw_coef' in directory:
	shutil.rmtree(directory)

for i in mat_coef_grid:
    for j in rmsd:
	directory = 'Matw_coef_'+str(i)+'_copy/rmsd'+str(j)
	os.system('mkdir -p '+directory)
	os.system('cp MR_pip.py'+' '+refl_file_input+' '+pdb_file_input+' '+seq_file_input+' '+directory)
	os.chdir("./"+directory)
    	os.system('bsub -q psanaq -n 12 -o %J.log python MR_pip.py -rfl '+refl_file_input+' -pdbE1 '+pdb_file_input+' -seq1 '+seq_file_input+' -idenE1 '+str(j)+' -errtE1 rmsd -c '+str(i)+' -res '+args.resolution+' -labin '+data_labels)
	os.chdir("../..")

