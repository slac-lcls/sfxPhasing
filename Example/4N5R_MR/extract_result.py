import pickle
import sys
import subprocess
import os
import argparse
import re
import json
import ast
import pandas as pd
import csv
#os.chdir('./GPCR_result')

unhidden_dir = []
for i in os.listdir(os.getcwd()):

    if len(i.split('.'))==1:
        unhidden_dir.append(i)

        
rmsdscann = []
files = []
for i in unhidden_dir:
    os.chdir(i)
    for j in os.listdir(os.getcwd()):
        rmsdscann.append(j)
        os.chdir(j)
        files.append(os.listdir(os.getcwd()))
        
        os.chdir('../')
    os.chdir('../')
    
result = []
for i in unhidden_dir:
    os.chdir(i)
    #print (unhidden_dir[i])
    for j in ['rmsd0.5','rmsd0.6','rmsd0.7','rmsd0.8','rmsd0.9','rmsd1.0','rmsd1.1','rmsd1.2','rmsd1.3','rmsd1.4','rmsd1.5','rmsd1.6','rmsd1.7','rmsd1.8','rmsd1.9','rmsd2.0']:
        os.chdir(j) #rmsdlist
        #print(rmsdscann[j])
        for k in os.listdir('.'):

            if k.split('.')[-1] == 'log':
                #print(k)
                #print('*')
                content=open(k, "r").read()
                component = 0
                result.append(i)
                result.append(j)
                for g in content.split('\n'):
                    if 'component TFZ' in g:
                        result.append(g)
                        component += 1
                    
                    if 'Run time' in g:
                        result.append(g)

                result.append('copy_got:'+str(component))
                result.append('####')
                component = 0

        os.chdir('..')
    os.chdir('..')

for i in range(len(result)):
    result[i] = result[i].replace(" ","")
    
list_out = []
list_in = []
result_final = []
for i in range(len(result)):
    list_in.append(result[i])
    if result[i] == '####':
        list_out.append(list_in)
        list_in = []
        
result_final = list_out

request_copy = []
rmsd = []
run_time= []
copy_got = []
TFZ = []
for i in result_final:

    one_set_TFZ = []
    for j in i:
        if j.split('.')[0].isdigit():
            one_set_TFZ.append(j.split('=')[1])
    TFZ.append(one_set_TFZ)       
    copy_got.append(i[-2].replace('copy_got:',''))
    run_time.append(i[-3].replace('Runtime:',''))
    rmsd.append(i[1].replace('rmsd',''))
    request_copy.append(i[0].replace('4RW2','').replace('copy',''))
    
TFZ_result=[]
string = '/'
for i in TFZ:
    for j in range(len(i)):
        string =string+ i[j]+'/'
        
    TFZ_result.append(string)
    string = '/'
    
Success=[]
for i in range(len(TFZ_result)):
    check = 0
    copies = 0
    for j in TFZ_result[i].split('/'):
        if j.replace('.','',1).isdigit():
            copies+=1
            if float(j) > 8:
                check+=1
                
    if check == copies:
        Success.append('Y')

    else:
        Success.append('N')
        

result_dictionary={'Requested Copy':request_copy, 'RMSD':rmsd, 
                   'TFZ':TFZ_result, 'Success':Success, 'Run Time':run_time, 'Copy Got':copy_got}


df = pd.DataFrame(data=result_dictionary)

df.to_csv(r'4RW2_result.csv')

