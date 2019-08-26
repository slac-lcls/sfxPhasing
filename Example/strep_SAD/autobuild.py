from __future__ import print_function
import sys
import subprocess
import os
import argparse
import re
import json
import ast

if os.path.isfile("autobuild.eff"):
    os.remove("autobuild.eff")

parser= argparse.ArgumentParser()

parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-seq","--sequence-file",help = "input the sequence file", type = str)
parser.add_argument("-rfff", "--r-free-flag-fraction-parameter", default = 0.1, help='enter the r free flag fraction', type = float)
parser.add_argument("-nproc", "--number-of-processors", default = 1, help='enter the number of processors', type = int)


args = parser.parse_args()

if args.reflection_mtz:
    reflectionFile = args.reflection_mtz
else:
    print('Please pass reflection file')

if args.sequence_file:
    sequenceFile = args.sequence_file
else:
    print('Please pass sequence file')


if args.r_free_flag_fraction_parameter:
    rFreeFlagFraction = args.r_free_flag_fraction_parameter

if args.number_of_processors:
    numberOfProcessors = args.number_of_processors

os.system('phenix.autobuild --show_defaults > my_autobuild.eff')

list = []
with open('my_autobuild.eff','r') as f:
    for line in f:
        list.append(line.rstrip('\n'))

list[0] = '#'
################################## Modify mtzfile 
print('#!/bin/csh -f',file=open("mtz_label_modification.sh", "a"))
print("mtzutils hklin2 strep.mtz hklin1 "+reflectionFile+" hklout test.mtz << eof",file=open("mtz_label_modification.sh", "a"))
print("merge",file=open("mtz_label_modification.sh", "a"))
print("eof",file=open("mtz_label_modification.sh", "a"))
print("#",file=open("mtz_label_modification.sh", "a"))
os.system("sh mtz_label_modification.sh")

reflectionFile = 'test.mtz'
########################################## Extract MTZ Information#################################
process = subprocess.Popen('phenix.mtz.dump '+reflectionFile, 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()

split_out=out.splitlines()

mylist = []
for i in range(len(split_out)):
    item = split_out[i].decode("utf-8")
    if item != '':
        mylist.append(split_out[i].decode("utf-8"))

for i in range(len(mylist)):
    if 'Space group from matrices:' in mylist[i]:
        SPACE_GROUP = mylist[i].replace('(',':').split(':')[1]

    elif 'Unit cell' in mylist[i]:
        CELL = mylist[i].replace('(',':').replace(')',':').split(':')[2].replace(',','')
#############################################################################################################

for i in range(len(list)):
    if ' data ' in list[i]:
        original = list[i].split('=')[1]
        list[i] = list[i].replace(original, ' "'+reflectionFile+'"')

    elif ' seq_file ' in list[i]:
        original = list[i].split('=')[1]
        list[i] = list[i].replace(original, ' "'+sequenceFile+'"')

    elif 'unit_cell' in list[i]:
        original = list[i].split('=')[1]
        list[i] = list[i].replace(original, ' '+CELL)

    elif 'space_group' in list[i]:
        original = list[i].split('=')[1]
        list[i] = list[i].replace(original, SPACE_GROUP)
    

    elif 'nproc' in list[i]:
        original = list[i].split('=')[1]
        list[i] = list[i].replace(original, ' '+str(numberOfProcessors))

    elif 'r_free_flags_fraction' in list[i]:
        original = list[i].split('=')[1]
        list[i] = list[i].replace(original, ' '+str(rFreeFlagFraction))

    elif 'clean_up' in list[i]:
        original = list[i].split('=')[1]
        list[i] = list[i].replace(original, ' True ')

for my_string in list:
    print(my_string,file=open("autobuild.eff", "a"))

os.system('phenix.autobuild autobuild.eff')