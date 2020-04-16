#!/usr/bin/env python
# Imports
import os, sys, re, subprocess, tqdm
import numpy as np
import pandas as pd
import datetime
import glob

# Example input
#python scrape.py 14601 72_RL_101_200 > creates 72_RL_101_200_date/time.pkl
 
# Arguments
project = sys.argv[1]                                             # Project Number
description = sys.argv[2]                                         # Description -i.e. Ligand_Set_L/RL_start_end
whole_dataset = pd.read_pickle(f"{re.split(r'[_-]', f'{description}')[0]}.pkl")
whole_project = whole_dataset.loc[whole_dataset['project'] == int(project)]  # extract info for project being scraped

try:                                                              # This loads in the most recently scraped df if it exists 
    previous_df = pd.read_pickle(max(glob.glob(f'scraped_data/{description}*') , key=os.path.getctime))
    chk_df = previous_df.filter(['run','clone','gen'])  # This will pull out run, clone, gen from old df
except Exception as e:                                            # if there is no old df...
    previous_df,chk_df = None, None
columns = ['description','identity','project','run','clone','gen','wl_increment'] + [f'lam_{x}' for x in range(39)] + ['free_energy'] # Columns for df

if not os.path.exists('scraped_data'):
    os.makedirs('scraped_data')

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

# server-dependent variables
hostname_paths = {'vav3.ocis.temple.edu':'/home/server/server2/data/SVR166219',
                  'vav4.ocis.temple.edu':'/home/server/server2/data/SVR166220',
                  'vav15.ocis.temple.edu':'/data/SVR2616698069',
                  'vav16.ocis.temple.edu':'/data/SVR2616698070',
                  'folder-prod-001':'/data/SVR1163805190', #avast1
                  'folder-prod-002':'/data/SVR1163805191', #avast2
                  'fah5':'/data/SVR679057516',
                  'aws2':'/data/SVR51748107'}
hostname = subprocess.check_output('hostname', shell=True).decode().strip('\n')
path = f'{hostname_paths[hostname]}/PROJ{project}'
    
# scrape the new data
data = []
for run in tqdm.tqdm(range(runs)):
    run_info = whole_project.loc[whole_project['run'] == run]
    for clone in range(clones):
        for gen in range(gens):
            if not os.path.exists(f'{path}/RUN{run}/CLONE{clone}/results{gen}/md.log'):    # This breaks the loop if the gen does not exist
                break
            if chk_df is not None:                                                       # This checks if run,clone,gen exists in old df
                if ((chk_df['run'] == run) & (chk_df['clone'] == clone) & (chk_df['gen'] == gen)).any() == True:
                    continue
            else:
                pass
            # this averages the last 100 entries in md.log
            try:
                cmd = f'tail -n 6208 {path}/RUN{run}/CLONE{clone}/results{gen}/md.log | grep -B 39 "40  0.000  1.000  1.000"'
                log = [line.split() for line in subprocess.check_output(cmd, shell=True).decode().split('\n')[:-1]]
                free_energies = [np.average([float(line[5]) for line in log if line[0] == str(lam)]) for lam in range(1,41)]
            except Exception as e:
                print(f'Bad md.log file! {path}/RUN{run}/CLONE{clone}/results{gen}/md.log. Skipping to next gen, I guess?')

            # SMILES...yr on Camera ( writes SMILES string or identity)  
            identity = run_info['identity'].values[0]

            try: # get WL-increment at end of md.log
                cmd = f'tail -n 1000 {path}/RUN{run}/CLONE{clone}/results{gen}/md.log | grep "increment" | tail -n 1'
                wl_increment = float(subprocess.check_output(cmd, shell=True).decode().split()[3])
                data.append([description, identity, project, run, clone, gen, wl_increment] + free_energies)
            except Exception as e:
                print(f'Bad md.log file2! {path}/RUN{run}/CLONE{clone}/results{gen}/md.log. Skipping to next gen, I guess?')
                continue

if chk_df is not None:
    df = previous_df.append(pd.DataFrame(data, columns=columns), ignore_index=True)     # This appends old df with new info
    df = df.sort_values(['run', 'clone', 'gen'], ascending=[True, True, True])          # This sorts run,clone,gen in sequential order
    df = df.reset_index(drop=True)                                                      # This resets the index to be sequential
else:                                                                                   
    df = pd.DataFrame(data, columns=columns)                                            # If no old df, make new one 
datetime_object = datetime.datetime.now()
ts = datetime_object.strftime("%d%b%Y_%H-%M-%S")
df.to_pickle(f'scraped_data/{description}_{ts}.pkl')

