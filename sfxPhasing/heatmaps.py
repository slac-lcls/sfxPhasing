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
import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#============#
# MR PHASING #
#============#

def Heat_map_prob(file,copy, upper_bound, lower_bound):
    raw_results = []
    with open(file,'r') as f:
        for line in f:
            raw_results.append(line.rstrip('\n'))

    result = []
    for line in raw_results:
        if "/Request_"+str(copy)+"_copy/" in line:
            result.append(line)

    resolution = []
    rmsd = []
    p_total = []
    for line in result:
        rmsd.append(float(line.split('rmsd')[-1].split('/')[0]))
        resolution.append(float(line.split('resolution')[-1].split(' ')[0]))
        p_total.append(float(line.split('P(total)=')[-1].split(' ')[0]))

    result_summary = pd.DataFrame(
        {'Total Probability': p_total,
         'rmsd': rmsd,
         'resolution': resolution
        })

    result_summary=result_summary.pivot(index='resolution', columns='rmsd', values='Total Probability')
    result_summary.to_latex(index=True)

    # Change Paramters of plot

    #sns.set(font_scale=1.8)
    plt.figure(figsize=(15, 10),dpi=400)
    sns.set(font_scale=1.8)
    sns_plot = sns.heatmap(result_summary, vmax=str(upper_bound),vmin=str(lower_bound),square=True, annot=False,
                 linewidths=.5 , xticklabels = 2,yticklabels = 2, cmap = sns.cm.rocket_r)
    sns_plot.set_xticklabels(sns_plot.get_xmajorticklabels(), fontsize = 20)
    sns_plot.set_yticklabels(sns_plot.get_ymajorticklabels(), fontsize = 20)

    #ax = sns.heatmap(result_summary, xticklabels=rmsd, yticklabels=resolution)
    #ax.set_xticks(ax.get_xticks()[::3])
    #ax.set_xticklabels(xlabels[::3])
    #ax.set_yticks(ax.get_yticks()[::3])
   # ax.set_yticklabels(ylabels[::3])

    plt.xlabel('R.m.s.d'+'$(\AA)$', fontsize = 25)
    plt.ylabel('Resolution'+'$(\AA)$', fontsize = 25)
    #plt.locator_params(axis='y', nbins=10)
    #plt.locator_params(axis='x', nbins=10)

    plt.savefig('GPCR_'+str(copy)+'_copy.eps', dpi = 300, metadata = 'eps',bbox_inches = 'tight',pad_inches = 0.1)
    plt.show()
    #plt.title('GPCR', fontsize = 20) # title with fontsize 20
    #plt.figure(figsize=(15, 15))





def Heat_map_TFZ(file,copy, upper_bound, lower_bound):
    raw_results = []
    with open(file,'r') as f:
        for line in f:
            raw_results.append(line.rstrip('\n'))

    result = []
    for line in raw_results:
        if "/Request_"+str(copy)+"_copy/" in line:
            result.append(line)

    resolution = []
    rmsd = []
    TFZ = []
    for line in result:
        rmsd.append(float(line.split('rmsd')[-1].split('/')[0]))
        resolution.append(float(line.split('resolution')[-1].split(' ')[0]))
        TFZ.append(float(line.split('TFZ = ')[-1].split(' ')[0]))

    result_summary = pd.DataFrame(
        {'TFZ': TFZ,
         'rmsd': rmsd,
         'resolution': resolution
        })

    result_summary=result_summary.pivot(index='resolution', columns='rmsd', values='TFZ')

    # Change Paramters of plot
    plt.figure(figsize=(15, 10),dpi=400)
    sns.heatmap(result_summary, vmax=str(upper_bound),vmin=str(lower_bound),square=True, annot=False,
                linewidths=.5 , xticklabels = 2,yticklabels = 2, cmap = sns.cm.rocket_r)
    plt.xlabel('R.m.s.d'+'$(\AA)$', fontsize = 20)
    plt.ylabel('Resolution'+'$(\AA)$', fontsize = 20)
    plt.savefig('lyso_'+str(copy)+'_copy.eps', dpi = 300, metadata = 'eps')
    plt.show()



#=============#
# SAD PHASING #
#=============#

def Heat_map_prob_SAD(FILE,THRESHOLD, upper_bound, lower_bound, DSUL ):
    raw_results = []
    with open(FILE,'r') as f:
        for line in f:
            raw_results.append(line.rstrip('\n'))

    result_dsul = []
    result_thre = []
    if DSUL != 0:
        for line in raw_results:
            if "DSUL"+str(DSUL) in line:
                result_dsul.append(line)
    else:
        for line in raw_results:
            result_dsul.append(line)

    for line in result_dsul:
        if "threshold"+str(THRESHOLD) in line:
            result_thre.append(line)

    resolution = []
    atom_number = []
    R_free = []
    for line in result_thre:
        resolution.append(float(line.split('resolution')[-1].split('/')[0]))
        atom_number.append(float(line.split('atom_number')[-1].split('/')[0]))
        if line.split('R_free:')[-1].split('/')[0] == ' ':
            R_free.append(round(float(line.split('R_free:')[-1].split(' ')[0]),3))
        else:
            R_free.append(round(float(line.split('R_free:')[-1].split('/')[0]),3))
    result_summary = pd.DataFrame(
        {'R_free': R_free,
         'atom_number': atom_number,
         'resolution': resolution
        })

    result_summary=result_summary.pivot(index='resolution', columns='atom_number', values='R_free')
    result_summary.to_latex(index=True)

    # Change Paramters of plot

    #sns.set(font_scale=1.8)
    plt.figure(figsize=(15, 10),dpi=400)
    sns.set(font_scale=1.8)
    sns_plot = sns.heatmap(result_summary, vmax=str(upper_bound),vmin=str(lower_bound),square=True, annot=True,
                 linewidths=.5 , xticklabels = 2,yticklabels = 2,fmt='.3g',annot_kws={"size": 15})
    sns_plot.set_xticklabels(sns_plot.get_xmajorticklabels(), fontsize = 20)
    sns_plot.set_yticklabels(sns_plot.get_ymajorticklabels(), fontsize = 20)

    #ax = sns.heatmap(result_summary, xticklabels=rmsd, yticklabels=resolution)
    #ax.set_xticks(ax.get_xticks()[::3])
    #ax.set_xticklabels(xlabels[::3])
    #ax.set_yticks(ax.get_yticks()[::3])
   # ax.set_yticklabels(ylabels[::3])

    plt.xlabel('Atom Number', fontsize = 25)
    plt.ylabel('Resolution'+'$(\AA)$', fontsize = 25)
    #plt.locator_params(axis='y', nbins=10)
    #plt.locator_params(axis='x', nbins=10)

    #plt.savefig('GPCR_'+str(copy)+'_copy.eps', dpi = 300, metadata = 'eps',bbox_inches = 'tight',pad_inches = 0.1)
    plt.show()
    #plt.title('GPCR', fontsize = 20) # title with fontsize 20
    #plt.figure(figsize=(15, 15))



