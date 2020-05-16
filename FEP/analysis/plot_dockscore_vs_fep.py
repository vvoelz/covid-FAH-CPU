#!/usr/bin/env python

import pickle
import pandas as pd
import matplotlib.pyplot as plt

db1 = pd.read_pickle('MS0326_results_1-5660.pkl')

fig, ax = plt.subplots()
ax.errorbar(db1['score'], db1['febkT']*2.49, yerr=db1['error']*2.49, fmt='.', ecolor='red')
ax.set_xlabel('Docking Score (Hybrid2)', fontsize=16)
ax.set_ylabel('FEP affinity (KJ/mol)', fontsize=16)
fig.savefig('docking_score_vs_fep.pdf')
