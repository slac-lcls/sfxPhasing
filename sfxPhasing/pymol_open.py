import sys, time, os
import argparse
import numpy as np
import pymol
from pymol import cmd

parser = argparse.ArgumentParser()
parser.add_argument("-map", "--reflection-mtz",help = "input the mtz file of reflection", type = str)
parser.add_argument("-pdb", "--pdb-file", help = "inpur the pdb file of reflection", type = str)
parser.add_argument("-r", "--radius", help = "input the radius you want to check the eletron density", type = float)

args = parser.parse_args()

if args.reflection_mtz:
    reflectionFile = args.reflection_mtz
    map_name = args.reflection_mtz.replace(".mtz","")
else:
    print("Please input the mtz file")
    sys.exit()
    
if args.pdb_file:
    pdbFile = args.pdb_file
    pdb_name = args.reflection_mtz.replace(".pdb","")
else:
    print('Please input the pdb file')
    sys.exit()
    
if args.radius:
    radius_check = args.radius
else:
    print("Input the radius you want to check")
    sys.exit()
    

pymol.finish_launching(["pymol","-q"])

last_m = '0'
last_p = '0'

while(1 == 1):
    #reflectionFile = 'strep_009.mtz'
    #pdbFile = 'strep_009.pdb'
    current_m = os.stat(reflectionFile)[8]
    current_p = os.stat(pdbFile)[8]
    
    if current_m != last_m or current_p != last_p:
        last_m = current_m
        last_p = current_p
        try:
            cmd.delete('my_map.2fofc')
            cmd.delete('my_map.fofc')
            cmd.load(pdbFile, 'overall_best')
            cmd.load_mtz(reflectionFile, prefix = 'my_map')#, amplitudes = 'unknown/unknown291019/REFM_FWT', 
                        #phases = 'unknown/unknown291019/REFM_PHWT',reso_low = 23.96, reso_high = 1.9)
            #cmd.map_double("my_map")
            cmd.map_double("my_map.2fofc")
            print("upload finished")
            cmd.hide("everything")
            cmd.show("sticks")
            cmd.set_bond('stick_radius','0.1', 'overall_best')
            
            cmd.center()
            
            cmd.select('site', 'br. all within '+str(radius_check)+' of center')
            
            cmd.isomesh('map','my_map.2fofc', 1.0, 'site', carve = 1.6)
            
            cmd.color('gray40','map')
            
            cmd.set('mesh_width','0.1')
                
            cmd.bg_color('white')
            
            cmd.set("ray_trace_fog","0")
            cmd.set("depth_cue","0")
            cmd.set("ray_shadows","off")
            cmd.ray(1024,1024)
            
            cmd.show("mesh","map")

            cmd.zoom('site')
            
        except:
            time.sleep(10)
            
    else:
        pass
            
            
            
            
            
            
