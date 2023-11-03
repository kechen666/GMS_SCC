import pandas as pd
#import subprocess
import os
import multiprocessing as mp
import numpy as np
from sys import argv

header = ['edge']
results_values = ['stateA', 'stateB', 'CGI_dG', 'CGI_err', 'BAR_dG', 'BAR_err', 'Forward_KS', 'Reverse_KS']
for leg in ['protein', 'water']:
    for value in results_values:
        header.append('_'.join([leg, value]))    
#print(header)


#with open('good_list_02.txt', 'r') as calculated_file:
#    calculated_edges = [i.strip('\n').strip(' ') for i in calculated_file.readlines()]

def parse_results(results_file):
    with open(results_file, 'r') as rfile:
        rfile_lines = rfile.readlines()
        #print(rfile_lines)
        for i,line in enumerate(rfile_lines):
            if 'Number of forward (0->1) trajectories' in line:
                num_forward = int(line.split()[-1])
            if 'Number of reverse (1->0) trajectories' in line:
                num_reverse = int(line.split()[-1])
            if 'CGI: dG =' in line:
                cgi_dg = float(line.split()[-2])
            if 'CGI: Std Err (bootstrap) =' in line:
                cgi_err = float(line.split()[-2])
            if 'Forward: gaussian quality =' in line:
                forward_ks = rfile_lines[i+1].split()[2]
            if 'Reverse: gaussian quality =' in line:
                reverse_ks = rfile_lines[i+1].split()[2]
            if 'BAR: dG =' in line:
                bar_dg = float(line.split()[-2])
            if 'BAR: Std Err (bootstrap)' in line:
                bar_err = float(line.split()[-2])
    return([num_forward, num_reverse, cgi_dg, cgi_err, bar_dg, bar_err, forward_ks, reverse_ks])

def calc_results(edge_list, out_results, calculate = True):
    for edge in edge_list:
        results_list = [edge]
        for leg in ['protein', 'water']:
            if calculate == True:
                print('Calculating values from '+leg)
                os.system(' '.join(['pmx', 'analyse',\
                    '-fA', '/'.join([edge,leg,'stateA/run?/nes_??.xvg']),\
                    '-fB', '/'.join([edge,leg,'stateB/run?/nes_??.xvg']),\
                    '-oA', '/'.join([edge,leg,'integA.dat']),\
                    '-oB', '/'.join([edge,leg,'integB.dat']),\
                    '-o',  '/'.join([edge,leg,'results.txt']),\
                    '-n', '100', '-t', '298.15', '-w', '/'.join([edge,leg,'wplot.png']),\
                    '--units', 'kcal', '--quiet']))
            results_list = results_list + parse_results('/'.join([edge,leg,'results.txt']))
        out_results.append(results_list)
    

edge = argv[1]
out_results = []

calc_results([edge],out_results, calculate = False)

results_df = pd.DataFrame(out_results, columns=header)
results_df['calculated_ddG_BAR'] = results_df['protein_BAR_dG'] - results_df['water_BAR_dG']
results_df['calculated_ddG_CGI'] = results_df['protein_CGI_dG'] - results_df['water_CGI_dG']
print('Your ΔΔG of binding in the '+edge+' system is '+ str(round(float(results_df['calculated_ddG_CGI'].iloc[0]), 2))+' kcal/mol according to the CGI method.')
print('Your ΔΔG of binding in the '+edge+' system is '+ str(round(float(results_df['calculated_ddG_BAR'].iloc[0]), 2))+' kcal/mol according to the BAR method.')
for state in ['protein_stateA', 'protein_stateB', 'water_stateA', 'water_stateB']:
    if results_df[state].iloc[0] != 256:
        print('WARNING! Not all runs completed successfully!')
        print(state+ ' has only '+ str(results_df[state].iloc[0])+' completed runs')
    else:
        print('You have sucessfully completed all runs in '+ state)


#print(results_df[['calculated_ddG', 'protein_BAR_dG', 'protein_BAR_err', 'water_BAR_dG', 'water_BAR_err',  'protein_stateA', 'protein_stateB', 'water_stateA', 'water_stateB']])


