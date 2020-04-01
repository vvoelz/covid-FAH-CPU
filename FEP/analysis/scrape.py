#!/usr/bin/env python

import os, sys, re, subprocess, tqdm
import numpy as np
import pandas as pd
import datetime

project = sys.argv[1]
description = sys.argv[2]
columns = ['description','SMILES','project','run','clone','gen','wl_increment'] + [f'lam_{x}' for x in range(39)] + ['free_energy']

# pull some info out of project.xml
with open(f'/home/server/server2/projects/Gromacs/p{project}/project.xml') as f:
    lines = [line.strip('\n').split() for line in f.readlines()]
for line in lines:
    try:
        if 'runs' in line[0]:
            runs = int(re.findall("\d+", line[1])[0])
        elif 'clones' in line[0]:
            clones = int(re.findall("\d+", line[1])[0])
        elif 'gens' in line[0]:
            gens = int(re.findall("\d+", line[1])[0])
        elif 'simulation' in line[1]:
            gen_time = float(line[4])
    except Exception as e:
        continue

path = f'../SVR2616698070/PROJ{project}'
if sys.argv[3] != None:
    previous_df = pd.read_pickle(sys.argv[3])
    chk_df = previous_df.filter(['run','clone','gen'])
else:
    pass

data = []
for run in tqdm.tqdm(range(runs)):
    for clone in range(clones):
        for gen in range(gens):
            if not os.path.exists(f'{path}/RUN{run}/CLONE{clone}/results{gen}/md.log'):
                break
            elif sys.argv != None: 
                if ((chk_df['run'] == run) & (chk_df['clone'] == clone) & (chk_df['gen'] == gen)).any() == True:
                    continue
            else:
                pass
# this averages the last 100 entries in md.log
            cmd = f'tail -n 6208 {path}/RUN{run}/CLONE{clone}/results{gen}/md.log | grep -B 39 "40  0.000  1.000  1.000"'
            log = [line.split() for line in subprocess.check_output(cmd, shell=True).decode().split('\n')[:-1]]
            free_energies = [np.average([float(line[5]) for line in log if line[0] == str(lam)]) for lam in range(1,41)]

#            cmd = f'tail -n 6208 {path}/RUN{run}/CLONE{clone}/results{gen}/md.log | grep "increment"'
#            log = [line.split() for line in subprocess.check_output(cmd, shell=True).decode().split('\n')[:-1]]
#            wc_increment = np.average([float(line[3]) for line in log])

# SMILES...yr on Camera ( writes SMILES string)  
            cmd = f'head -n1 ../projects/p{project}/RUN{run}/*_h.sdf'
            SMILES = subprocess.check_output(cmd,shell=True)

# this uses just the last entry of md.log
#            cmd = f'tail -n 1000 {path}/RUN{run}/CLONE{clone}/results{gen}/md.log | grep -B 40 "40  0.000  1.000  1.000" | tail -n 40'
#            free_energies = [line.split() for line in subprocess.check_output(cmd, shell=True).decode().split('\n')[:-1]]

            cmd = f'tail -n 1000 {path}/RUN{run}/CLONE{clone}/results{gen}/md.log | grep "increment" | tail -n 1'
            wl_increment = float(subprocess.check_output(cmd, shell=True).decode().split()[3])
            data.append([description, SMILES, project, run, clone, gen, wl_increment] + free_energies)

if sys.argv[3] != None:
    df = previous_df.append(pd.DataFrame(data, columns=columns), ignore_index=True)
    df = df.sort_values(['run', 'clone', 'gen'], ascending=[True, True, True])
    df = df.reset_index()
else:
    df = pd.DataFrame(data, columns=columns)
datetime_object = datetime.datetime.now()
ts = datetime_object.strftime("%d%b%Y_%H-%M-%S")
df.to_pickle(f'avg_{description}_{ts}.pkl')

