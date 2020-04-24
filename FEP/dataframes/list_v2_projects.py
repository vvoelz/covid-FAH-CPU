#! /usr/bin/env python

import os, sys
import pandas as pd



usage = """Usage:

    python list_v2_projects.py [seriesname]

Will list the directories to convert to --> from  for v2 projects, where seriesname is e.g.  'MS0323_RL'
    
    Example: $$ python list_v2_projects.py MS0323_RL
"""

if len(sys.argv) < 2:
    print(usage)
    sys.exit(1)

seriesname = sys.argv[1]

dataframe_dir = '/home/server/git/covid-FAH-CPU/FEP/dataframes'

master_pkl_file = os.path.join(dataframe_dir, 'master_FEP.pkl')
df = pd.read_pickle(master_pkl_file)


# if you want to iterate over rows that match your selection do:
for index, row in df.loc[df['dataset'] == seriesname].iterrows():
    print(f"p{row['v1_project']}/RUN{row['v1_run']} --> p{row['v2_project']}/RUN{row['v2_run']}")

# you can pair multiple selections together in a similar loop:
#for index, row in df.loc[(df['dataset'] == '72_RL') & (df['v2_run'] < 100)].iterrows():
#    print(row)


