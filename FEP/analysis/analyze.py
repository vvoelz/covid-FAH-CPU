#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import tqdm, sys, glob, datetime, os
import re, subprocess

# in columns = [description  project  run  clone  gen  wl_increment  lam_0 ... lam38 free energy]
# out columns [run  ns_RL  ns_L  wl_increment_RL  wl_increment_L  free_energy  L39..L0  RL0..R39]
# Sample input: Ligset_start#_end# < this is the description variable

description = sys.argv[1] # should be like MS0326_1-1500 or AGG_10000-12000
project = sys.argv[2] # project number assoicated with RL-system matching description
lambdas = 40 # 20 lambdas to (de)couple each of coulomb/vdw == 40 total

# constraints
RL_gen_cutoff = 10 # only process clones that have more than N ns of sampling
L_gen_cutoff = 5
wl_increment_cutoff = 0.25 # only plot ligands with avg_RL_WL_increment < N
clone_cutoff = 2 # only include runs that have >=N clones that satisfy restraints (use 1-3)

# get FAH data directory from hostname
hostname = subprocess.check_output('hostname', shell=True).decode().strip('\n')
hostname_paths = {'vav3.ocis.temple.edu':'/home/server/server2/data/SVR166219',
                  'vav4.ocis.temple.edu':'/home/server/server2/data/SVR166220',
                  'vav15.ocis.temple.edu':'/data/SVR2616698069',
                  'vav16.ocis.temple.edu':'/data/SVR2616698070',
                  'folder-prod-001':'/data/SVR1163805190', #avast1
                  'folder-prod-002':'/data/SVR1163805191', #avast2
                  'fah5':'/data/SVR679057516',
                  'aws2':'/data/SVR51748107'}
server_path = f'{hostname_paths[hostname]}/PROJ{project}'

try:
    desc = re.split(r'[_-]', f'{description}')
    data_RL_file = max(glob.glob(f'scraped_data/{desc[0]}_RL_{desc[-2]}*{desc[-1]}*.pkl'), key=os.path.getctime)
    data_L_file = max(glob.glob(f'scraped_data/{desc[0]}_L_{desc[-2]}*{desc[-1]}*.pkl'), key=os.path.getctime)
    data_RL = pd.read_pickle(data_RL_file)
    data_L = pd.read_pickle(data_L_file)
    whole_dataset = pd.read_pickle(f'{desc[0]}.pkl') # master dataframe for dataset (ask matt if missing)
    print(f'Using most recent data files:\n\tRL Data: {data_RL_file}\n\tL Data: {data_L_file}')
#    print(f'{data_RL}\n{data_L}')
except IndexError as e:
    print(f'No pickled dataset found for description: {description}')

# make sure output directories are present
for dir in ['plots','results']:
    if not os.path.exists(dir):
        os.makedirs(dir)

## scrape pull info
temperature = 300.0
kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
with open(f'/home/server/server2/projects/Gromacs/p{project}/RUN0/prod.mdp') as f:
    lines = [line.strip('\n').split() for line in f.readlines()]
for line in lines:
    try:
        if 'pull-coord1-init' in line[0]:
            equilibrium_distance = float(line[2])
        elif 'pull_coord1_k' in line[0]:
            kspring = float(line[2])
    except Exception as e:
        continue

# process runs that have data for both RL/L
runs = set(list(data_RL['run'].values) + list(data_L['run'].values))
data, errors, clone_energies = [],[],[]
for run in tqdm.tqdm(runs):
    run_RL = data_RL.loc[data_RL['run'] == run]
    run_L = data_L.loc[data_L['run'] == run]

    # process L systems for last gen of each clone
    energies_L, wl_increment_L, ns_L = [],[],[]
    for clone_L in set(run_L['clone'].values):
        try:
            clone_L_data = run_L.loc[run_L['clone'] == clone_L]
            gen_L = clone_L_data.loc[clone_L_data['gen'] == clone_L_data['gen'].max()]
            if clone_L_data['gen'].max() < L_gen_cutoff:
                continue # skip clones with fewer than N gens or HIGH WL-increment
            wl_increment_L.append(gen_L['wl_increment'].values[0])
            ns_L.append(10+int(gen_L['gen'].values[0])*10)
            raw_energies = gen_L[gen_L.columns[-lambdas:]].values[0]
            energies_L.append([float(x) for x in raw_energies]) #- raw_energies[0] for x in raw_energies])
            clone_energies.append([run, 'L', ns_L[-1], wl_increment_L[-1], list(reversed(energies_L[-1]))])
        except Exception as e:
            print(f'L Exception: {run}, {clone_L}, {e}')
            continue
    try:
        avg_WL_increment_L = np.average(wl_increment_L)
        avg_ns_L = np.average(ns_L)
        avg_L_energies = [np.average(np.asarray(energies_L)[:,lam]) for lam in range(lambdas)]
    except Exception as e:
        continue

    # process RL systems for last gen of each clone
    energies_RL, wl_increment_RL, ns_RL, pull_energies = [],[],[],[]
    for clone_RL in set(run_RL['clone'].values):
        try:
            clone_RL_data = run_RL.loc[run_RL['clone'] == clone_RL]
            gen_RL = clone_RL_data.loc[clone_RL_data['gen'] == clone_RL_data['gen'].max()]
            if clone_RL_data['gen'].max() < RL_gen_cutoff or clone_RL_data['wl_increment'].min() > wl_increment_cutoff:
                continue # skip clones with fewer than N gens or HIGH WL-increment
            wl_increment_RL.append(gen_RL['wl_increment'].values[0])
            ns_RL.append(int(gen_RL['gen'].values[0])+1)
            raw_energies = gen_RL[gen_RL.columns[-lambdas:]].values[0]
            energies_RL.append([float(x) for x in raw_energies])
            clone_energies.append([run, 'RL', ns_RL[-1], wl_increment_RL[-1], energies_RL[-1]])
            # look at harmonic restraint of last gen, for each clone that fits the constraints
            with open(f"{server_path}/RUN{run}/CLONE{clone_RL}/results{int(gen_RL['gen'].values[0])}/pullx.xvg") as xvg:
                lines = [ line.strip('\n').split() for line in xvg.readlines()][20:][::10]
            displacements = [np.sqrt(float(line[1])**2 + float(line[2])**2 + float(line[3])**2) for line in lines]
            pull_energies.append(displacements) # only look at displacements for v1
#            pull_energies.append([beta*(kspring/2.0)*(displacement - equilibrium_distance)**2 for displacement in displacements])

        except Exception as e:
            print(f'RL Exception: {run}, {clone_RL}, {e}')
            continue

    try:
        avg_WL_increment_RL = np.average(wl_increment_RL)
        avg_ns_RL = np.average(ns_RL)
        avg_RL_energies = [np.average(np.asarray(energies_RL)[:,lam]) for lam in range(lambdas)]
        avg_free_energy = avg_RL_energies[-1] - avg_L_energies[-1]
        data.append([run, avg_ns_RL, avg_ns_L, avg_WL_increment_RL, avg_WL_increment_L, avg_free_energy] + list(reversed(avg_L_energies)) + list(avg_RL_energies))
    except Exception as e:
        continue

columns = ['run', 'ns_RL', 'ns_L', 'wl_increment_RL', 'wl_increment_L', 'free_energy'] + [f'L{lam}' for lam in reversed(range(lambdas))] + [f'RL{lam}' for lam in range(lambdas)]
clone_energies_columns = ['run', 'type', 'length', 'wl_increment', 'energies']
results = pd.DataFrame(data, columns=columns)
clone_energies = pd.DataFrame(clone_energies, columns=clone_energies_columns)
whole_project = whole_dataset.loc[whole_dataset['project'] == int(project)]
np.save(f'results/pull_{description}.npy',pull_energies)

print(f'*** Results are based on RL sampling > {RL_gen_cutoff} ns, L sampling > {L_gen_cutoff} ns, and a RL WL-increment < {wl_increment_cutoff}:')
### Saving Plots and Dataframe
good_results = []
good_results_columns = ['dataset','fah','identity','receptor','score','febkT','error','ns_RL','ns_L','wl_RL']
for run in tqdm.tqdm(results['run'].values):
    try:
        run_info = whole_project.loc[whole_project['run'] == run]
        clones_RL = clone_energies.loc[(clone_energies['run'] == run) & (clone_energies['type'] == 'RL')]
        clones_L = clone_energies.loc[(clone_energies['run'] == run) & (clone_energies['type'] == 'L')]
        result = results.loc[results['run'] == run]
        summary = list(result[result.columns[:6]].values[0]) # run ns_RL, ns_L, wl_increment_RL, wl_increment_L
        energies = result[result.columns[6:]].values[0]
        energies = [x - energies[0] for x in energies] # set first point to 0kT reference

        if len(clones_RL) < clone_cutoff: # skip runs which have fewer RL clones than the clone_cutoff
            continue

        clone_combinations = []
        fig, ax = plt.subplots()
        for lindex, lrow in clones_L.iterrows():
            for rindex, rrow in clones_RL.iterrows():
                clone_energy = lrow['energies'] + rrow['energies']
                clone_energy = [x - clone_energy[0] for x in clone_energy]
                ax.plot(range(len(clone_energy)), clone_energy)
                clone_combinations.append(clone_energy)
        energy_errors = np.std(clone_combinations, axis=0)
        ax.errorbar(range(len(energies)), energies, yerr=energy_errors)
        ax.axhline(0,linestyle='--') # draw horizontal lines at reference
        ax.axhline(summary[-1],linestyle='--') # and at free energy
        if result['wl_increment_RL'].values[0] > 0.2: # color complex lambdas based on WL-increment
            wl_color = 'red'
        elif result['wl_increment_RL'].values[0] > 0.1:
            wl_color = 'yellow'
        else:
            wl_color = 'green'
        ax.axvspan(0, lambdas-1, alpha=0.1, color='blue') # color first half of lambdas blue
        ax.axvspan(lambdas-1, lambdas*2, alpha=0.1, color=wl_color) # color second half of lambdas (consider changing this based on some other metric.)
        if summary[-1] < 0: # change color of horizontal free energy strip depending on sign
            nrg_color = 'green'
        else:
            nrg_color = 'red'
        ax.axhspan(0,summary[-1], alpha=0.5, color=nrg_color)
        ax.set_title(f"{desc[0]}: p{run_info['project'].values[0]}:R{run_info['run'].values[0]}, {run_info['identity'].values[0]}, ΔG = {summary[-1]:.3f}±{energy_errors[-1]:.2f}kT")
        plt.xlim(0,80)
        plt.xlabel('Ligand Decoupling → Receptor-Ligand Coupling')
        plt.ylabel('Free Energy (kT)')
        date = datetime.datetime.now()
        ts = date.strftime('%d%b%Y')
        plt.savefig(f'plots/{description}_{run}_{ts}.png')
        plt.close()

        good_results.append([desc[0],f"{server_path}/PROJ{run_info['project'].values[0]}/RUN{run_info['run'].values[0]}",run_info['identity'].values[0],run_info['receptor'].values[0],run_info['score'].values[0], result['free_energy'].values[0], energy_errors[-1], clones_RL['length'].values, clones_L['length'].values, clones_RL['wl_increment'].values])

    except Exception as e:
        print(f'Exception in computing energies/errors and plotting: {e}')
        continue

good_results = pd.DataFrame(good_results, columns=good_results_columns)
good_results = good_results.sort_values(by=['febkT'])
print(good_results)
good_results.to_pickle(f'results/{description}_{ts}.pkl')

# stuff I took out of plotting
#    results.to_pickle(f'results_{description}.pkl')
#    clone_energies.to_pickle(f'clone_energies_{description}.pkl')
#        print(f"Plotting RUN{run}\nRL_WL_increment: {result['wl_increment_RL'].values[0]:.3f}\nΔG_unb = {result['free_energy'].values[0]:.3f}±{energy_errors[-1]:.3f}kT\n")
#        ax.set_title(f'{desc[0]}_{run}: RL > {RL_gen_cutoff} ns, L > {L_gen_cutoff} ns, and RL-WLI < {wl_increment_cutoff:.3f}') #_RUN{summary[0]}: ΔG = {summary[-1]:.2f} ± {summary_error[-1]:.2f}kT')
#        label = f'RL/L: RUN{run} {int(summary[1])}/{int(summary[2])}ns {summary[3]:.3f}/{summary[4]:.3f}WL, ΔG = {summary[-1]:.3f}±{energy_errors[-1]:.2f}kT'
#        ax.annotate(label, xy = (range(len(energies))[-1], energies[-1]), xytext = (60, 100),
#              textcoords = 'offset points', ha = 'right', va = 'bottom',
#              bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
#                  arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
