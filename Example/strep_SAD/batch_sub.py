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

original_path = os.getcwd()
parser= argparse.ArgumentParser()

parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-seq","--sequence-file",help = "input the sequence file", type = str)
args = parser.parse_args()

if args.reflection_mtz:
   reflectionFile = args.reflection_mtz
   job_name = args.reflection_mtz.replace('.mtz','')
else:
   print('Please input the mtz file')

if args.sequence_file:
   sequenceFile = os.path.abspath(args.sequence_file)
else:
   print('Please input sequence')
   sys.exit()

resolution_range = np.arange(2.0, 4.1, 0.1)
atom_find = np.arange(2,9,1)

for resolution in resolution_range:
    for number in atom_find:
        directory = 'resolution'+str(resolution)+'/atom_number'+str(number)
        os.system('mkdir -p '+directory)
        os.system('cp Se_SAD_automation.py SHELX_script.py crank2_script.py autobuild.py '+sequenceFile+' '+reflectionFile+' '+directory)
        os.chdir("./"+directory)

        automation_cl = 'python Se_SAD_automation.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)
        os.system('bsub -q psanaq -n 12 -o %J.log '+automation_cl)
        os.chdir(original_path)
