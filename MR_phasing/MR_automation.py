import sys
import shutil
import subprocess
import os
import argparse
import re
import numpy as np

##################################### Parse file and parameters##############################################
# parse reflection file (mtz), model file (pdb), protein sequence file (pdb)
# parse a wanted resolution value

parser= argparse.ArgumentParser()

refl_file_input = ''
pdb_file_input = ''
seq_file_input = ''
data_labels = ''

parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-resl","--resolution", default = "2.5", help="define the resolution", type = str)
parser.add_argument("-pdb","--pdb-file",help="input the pdb file",type = str)
parser.add_argument("-seq","--sequence-file",help="input the sequence file",type = str)
parser.add_argument("")

args = parser.parse_args()

if args.reflection_mtz:
    refl_file_input = args.reflection_mtz
else:
    print('The reflection file is missing')
    sys.exit()

if args.pdb_file:
     pdb_file_input= args.pdb_file
else:
    print('The pdb file is missing')
    sys.exit()

if args.sequence_file:
    seq_file_input = args.sequence_file
else:
    print('The sequence file is missing')
    sys.exit()
#################################################################################################################


##################################### Get request copy number range ############################################
# Create the number of reasonable request copies by using phenix.xtriage to read the reflection file and sequence file
# If the read Matthew Coefficient = m, the requested copy is [m-1, m+1] with interaval of 1.
# Data label is also read here. The idea is to open one file and extract as much information as possible
process = subprocess.Popen('phenix.xtriage '+refl_file_input+' '+pdb_file_input+' '+seq_file_input, 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()

split_out=out.splitlines()
my_list=[]
for i in range(len(split_out)):
    my_list.append(split_out[i].decode("utf-8"))

for line in my_list:
    if 'Best guess' in line:
        string = line
        matthew_coefficient = int(re.findall('\d+', string )[0])

    if "Data labels" in line:
        data_labels = line.split(':')[1].replace(' ','')
        
        
if matthew_coefficient == 1:
    mat_coef_grid = range(1,3)
else:
    mat_coef_grid = range(matthew_coefficient-1,matthew_coefficient+2)
###############################################################################################################

#####################################Create the RMSD Grid###############################################
# define this arbitrarily
rmsd = np.linspace(0.5,2.0,16)
########################################################################################################

#########################################Get the resolution range#################################################
# Get resolution by extracting the highest resolution bound, say res_h, from the input mtz file
# Create the grid of [res_h, res_h+1.5) with interval of 0.1

# Using phenix.mtz.dump to read the reflection file 
process = subprocess.Popen('phenix.mtz.dump '+reflectionFile, 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()

split_out=out.splitlines()

for line in split_out:
    if 'Resolution range' in line:
        resolution = round(float(line.split(' ')[-1]),1)

# Create the grid of resolution
resolution_range = np.arange(resolution, resolution+1.5, 0.1)
##################################################################################################################

##################################### Implement the Grid and submit jobs ########################################
for directory in os.listdir(os.getcwd()):
    if 'Request_' in directory:
        shutil.rmtree(directory)

#create each grid as a directory
for i in mat_coef_grid:
    for j in rmsd:
        directory = 'Request_'+str(i)+'_copy/rmsd'+str(j)
        os.system('mkdir -p '+directory)
        os.system('cp MR_pip.py'+' '+refl_file_input+' '+pdb_file_input+' '+seq_file_input+' '+directory)
        os.chdir("./"+directory)
        os.system('bsub -q psanaq -n 12 -o %J.log python MR_pip.py -rfl '+refl_file_input+' -pdbE1 '+pdb_file_input+' -seq1 '+seq_file_input+' -idenE1 '+str(j)+' -errtE1 rmsd -c '+str(i)+' -res '+args.resolution+' -labin '+data_labels)
        os.chdir("../..")
#################################################################################################################
