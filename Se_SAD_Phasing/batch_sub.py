# THis is the top file of the SAD pipeline. 
# It is used to carry out the batch submission
# This .py will do several automation and grid definition, and parse
# these information to SAD-automation file.
#
from __future__ import print_function
import sys
import subprocess
import os
import argparse
import re
import json
import ast
import numpy as np
import shutil

#This is the original directory which is the 'root' directory of your results
original_path = os.getcwd()
parser= argparse.ArgumentParser()

parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-seq","--sequence-file",help = "input the sequence file", type = str)
parser.add_argument("-SFAC","--atom-type", help = "input the name of atom of this SAD (case insensitive)" , type = str)
args = parser.parse_args()

#parse reflection file
if args.reflection_mtz:
    reflectionFile = args.reflection_mtz
    job_name = args.reflection_mtz.replace('.mtz','')
else:
    print('Please input the mtz file')
    sys.exit()

if args.sequence_file:
    sequenceFile = os.path.abspath(args.sequence_file)
else:
    print('Please input sequence')
    sys.exit()
    
if args.atom_type:
    atomType = (args.atom_type).upper()
else:
    print('Please input the scattering atoms')
    sys.exit

###################### Setup the Grid Range #########################

print("This program utilize two softwares, CCP4i2 SHELXC/D and Crank2. Further refinement can be done by initiating Autobuild")

# Define the minimum distance between atoms and symmetry units. Users can change the definition in the parameter.json file in the 
# same directory.
with open('parameter.json') as json_file:
    data = json.load(json_file)
    for p in data['MIND1']:
        if p == atomType:
            mind_atom = str(data['MIND1'][atomType])
            
    for p in data['MIND2']:
        if p == atomType:
            mind_symm = str(data['MIND2'][atomType])
            
    for p in data["Low Resolution CutOff"]:
        if p == atomType:
            low_resolution_cut = str(data['Low Resolution CutOff'][atomType])
            
#Define the SE, S, and DSUL number from the sequence file
sequence=[]
with open(sequenceFile,'r') as f:
    for line in f:
        sequence.append(line.rstrip('\n'))
        
for line in sequence:
    if '>' in line:
        sequence.remove(line)

protein = ''.join(sequence)

single_S_or_SE_number=0
double_sulfur_number=0

for i in protein:
    if i=='M':
        single_S_or_SE_number+=1 #For S-SAD, M = S; For Se_SAD, M=SE
    elif i=='C':
        double_sulfur_number+=1
        
max_DSUL = double_sulfur_number/2
max_S = single_S_or_SE_number+double_sulfur_number
max_SE = single_S_or_SE_number

#######################################################
process = subprocess.Popen('phenix.mtz.dump '+'a2a_ccp4ifw.mtz', 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()

split_out=out.splitlines()

for line in split_out:
    if 'Resolution range' in line:
        resolution = round(float(line.split(' ')[-1]),1)
        

DSUL_range = range(1,max_DSUL+1)
resolution_range = np.arange(resolution, resolution+1.5, 0.1) #2.4 4.1 0.1
#atom_find = np.arange(10,13,1)#10,21,1
if atomType == 'SE':
    atom_find = range(max_SE/2, max_SE+1)
elif atomType == 'S':
    atom_find = range(max_S/2, max_S+1)

    
    
thre_range = np.linspace(0.2,0.5,4)


#consider for DSUL parameter
if max_DSUL > 0 and atomType == 'S':
    for dsul in DSUL_range:
        for thre in thre_range:
            for resolution in resolution_range:
                for number in atom_find:
                    directory = 'DSUL'+str(dsul)+'/threshold'+str(thre)+'/resolution'+str(resolution)+'/atom_number'+str(number)
                    os.system('mkdir -p '+directory)
                    os.system('cp Se_SAD_automation.py SHELX_script.py crank2_script.py autobuild.py '+sequenceFile+' '+reflectionFile+' '+directory)
                    os.chdir("./"+directory)

                    automation_cl = 'python Se_SAD_automation.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.3 -thre '+str(thre)+' -DSUL '+str(dsul)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm+' -lresl '+low_resolution_cut
                    os.system('bsub -q psanaq -n 12 -o %J.log '+automation_cl)
                    os.chdir(original_path)

#Do not consider DSUl parameter
elif max_DSUL == 0 or atomeType != 'S':
    print ('0 disulfide bond found or the atom type you are looking for is not Sulfur. Do not consider the grid search for disulfied number')
    for thre in thre_range:
        for resolution in resolution_range:
            for number in atom_find:
                directory = 'threshold'+str(thre)+'/resolution'+str(resolution)+'/atom_number'+str(number)
                os.system('mkdir -p '+directory)
                os.system('cp Se_SAD_automation.py SHELX_script.py crank2_script.py autobuild.py '+sequenceFile+' '+reflectionFile+' '+directory)
                os.chdir("./"+directory)

                automation_cl = 'python Se_SAD_automation.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.5 -thre '+str(thre)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm
                os.system('bsub -q psanaq -n 12 -o %J.log '+automation_cl)
                os.chdir(original_path)