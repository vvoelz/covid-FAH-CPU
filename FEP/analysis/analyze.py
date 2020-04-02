#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import tqdm, sys, glob

# in columns = [description  project  run  clone  gen  wl_increment  lam_0 ... lam38 free energy]
# out columns [run  ns_RL  ns_L  wl_increment_RL  wl_increment_L  free_energy  L39..L0  RL0..R39]

max_gen = True
data = []
description = '72_L_101_200'
data_RL = pd.read_pickle(glob.glob('avg_72_RL_101_200_*.pkl')[0])
data_L = pd.read_pickle(glob.glob('avg_72_L_101_200_*.pkl')[0])
lambdas = 40
RL_gen_cutoff = 15 # only process clones that have more than N ns of sampling
L_gen_cutoff = 10
wl_increment_cutoff = 0.2 # only plot ligands with avg_RL_WL_increment < N
skip_positive_G = True
columns = ['run', 'ns_RL', 'ns_L', 'wl_increment_RL', 'wl_increment_L', 'free_energy'] + [f'L{lam}' for lam in reversed(range(lambdas))] + [f'RL{lam}' for lam in range(lambdas)]
runs = min([data_RL['run'].max(), data_L['run'].max()]) + 1 # only process runs that have data for both RL and L (this can be done better)
data, errors = [],[]


for run in tqdm.tqdm(range(runs)):
    run_RL = data_RL.loc[data_RL['run'] == run]
    run_L = data_L.loc[data_L['run'] == run]

    # process L systems
    energies_L, wl_increment_L, ns_L = [],[],[]
    for clone_L in set(run_L['clone'].values):
        try:
            clone_L_data = run_L.loc[run_L['clone'] == clone_L]
            gen_L = clone_L_data.loc[clone_L_data['gen'] == clone_L_data['gen'].max()]
            if clone_L_data['gen'].max() < L_gen_cutoff: # we skip clones with fewer than N gens
                continue
            wl_increment_L.append(gen_L['wl_increment'].values[0])
            ns_L.append(int(gen_L['gen'].values[0]))
            raw_energies = gen_L[gen_L.columns[-lambdas:]].values[0]
            energies_L.append([x - raw_energies[0] for x in raw_energies])
        except Exception as e:
            print(f'L Exception: {run}, {clone_L}, {e}')
            continue
    try:
        avg_WL_increment_L = np.average(wl_increment_L)
        std_WL_increment_L = np.std(wl_increment_L)
        avg_ns_L = np.average(ns_L)
        std_ns_L = np.std(ns_L)
        avg_L_energies = [np.average(np.asarray(energies_L)[:,lam]) for lam in range(lambdas)]
        std_L_energies = [np.std(np.asarray(energies_L)[:,lam]) for lam in range(lambdas)]
    except Exception as e:
        continue

    # process RL systems
    energies_RL, wl_increment_RL, ns_RL = [],[],[]
    for clone_RL in set(run_RL['clone'].values):
        try:
            clone_RL_data = run_RL.loc[run_RL['clone'] == clone_RL]
            gen_RL = clone_RL_data.loc[clone_RL_data['gen'] == clone_RL_data['gen'].max()]
            if clone_RL_data['gen'].max() < RL_gen_cutoff: # we skip clones with fewer than N gens
                continue
            wl_increment_RL.append(gen_RL['wl_increment'].values[0])
            ns_RL.append(int(gen_RL['gen'].values[0]))
            raw_energies = gen_RL[gen_RL.columns[-lambdas:]].values[0]
            energies_RL.append([x - raw_energies[0] for x in raw_energies])
        except Exception as e:
            print(f'RL Exception: {run}, {clone_RL}, {e}')
            continue

    try:
        avg_WL_increment_RL = np.average(wl_increment_RL)
        std_WL_increment_RL = np.std(wl_increment_RL)
        avg_ns_RL = np.average(ns_RL)
        std_ns_RL = np.std(ns_RL)
        avg_RL_energies = [np.average(np.asarray(energies_RL)[:,lam]) for lam in range(lambdas)]
        std_RL_energies = [np.std(np.asarray(energies_RL)[:,lam]) for lam in range(lambdas)]
        avg_free_energy = avg_RL_energies[-1] - avg_L_energies[-1]
        std_free_energy = np.sum(std_RL_energies + std_L_energies)

        data.append([run, avg_ns_RL, avg_ns_L, avg_WL_increment_RL, avg_WL_increment_L, avg_free_energy] + list(reversed(avg_L_energies)) + list(avg_RL_energies))
        errors.append([run, std_ns_RL, std_ns_L, std_WL_increment_RL, std_WL_increment_L, std_free_energy] + list(reversed(std_L_energies)) + list(std_RL_energies))

    except Exception as e:
        continue

results = pd.DataFrame(data, columns=columns)
errors = pd.DataFrame(errors, columns=columns)
results.to_pickle(f'{description}_results.pkl')
errors.to_pickle(f'{description}_errors.pkl')

good_results = results.loc[results['wl_increment_RL'] < wl_increment_cutoff]
good_results_errors = errors.loc[results['wl_increment_RL'] < wl_increment_cutoff]

if skip_positive_G:
    good_results_errors = good_results_errors.loc[good_results['free_energy'] < 0]
    good_results = good_results.loc[good_results['free_energy'] < 0]

print(f'*** Results based on RL sampling > {RL_gen_cutoff} ns, L sampling > {L_gen_cutoff} ns, and a RL WL-increment < {wl_increment_cutoff}:')
print(good_results)

### example plots
for run in good_results['run'].values:
    try:
        result = good_results.loc[good_results['run'] == run]
        result_error = good_results_errors.loc[good_results['run'] == run]
        print(f"Plotting RUN{run}\nRL_WL_increment: {result['wl_increment_RL'].values[0]:.3f}\nΔG_unb = {result['free_energy'].values[0]:.3f}kT\n")
        summary = list(result[result.columns[:6]].values[0]) # run ns_RL, ns_L, wl_increment_RL, wl_increment_L
        summary_error = list(result_error[result_error.columns[:6]].values[0])

        energies = result[result.columns[6:]].values[0]
        energies = [x - energies[0] for x in energies] # set first point to 0
        energy_errors = result_error[result_error.columns[6:]].values[0]
        energy_errors = [sum(energy_errors[0:x + 1]) for x in range(len(energy_errors))]

        plt.scatter(range(len(energies)), energies)
        plt.plot(range(len(energies)), energies)
        plt.errorbar(range(len(energies)), energies, yerr=energy_errors)
        plt.title(f'MS0323_RUN{summary[0]}: ΔG = {summary[-1]:.3f} ± {summary_error[-1]:.3f}kT')  
        plt.xlabel('Ligand Decoupling --> Receptor-Ligand Coupling')
        plt.ylabel('Free Energy (kT)')
        plt.savefig(f'FEP_{run}.png')
#        plt.close()
    except Exception as e:
        print(f'Exception in computing energies/errors and plotting: {e}')
        continue
