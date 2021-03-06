This directory contains dataframes for each existing dataset as well as a jupyter notebook that was used to construct the master dataframe,
"master_FEP.pkl". This contains information on the identity, score, and project/run for each v1/v2/v3 for each ligand of each dataset.
To access this information, use:

```
import pandas as pd
df = pd.read_pickle('master_FEP.pkl')
```

You can access/print a certain subset of this by any of the columns, so let's list columns with:
```
df.columns
```

Then you can look at all of the 'MS0323' dataset with:
```
df.loc[df.dataset.str.contains('MS0323')]
```

or just the RL systems of the '387' dataset with:
```
df.loc[df['dataset'] == '387_RL']
```

if you want to iterate over rows that match your selection do:
```
for index, row in df.loc[df['dataset'] == '72_RL'].iterrows():
    print(f"cp p{row['v1_project']}/RUN{row['v1_run']} p{row['v2_project']}/RUN{row['v2_run']}")
```

you can pair multiple selections together in a similar loop:

```
for index, row in df.loc[(df['dataset'] == '72_RL') & (df['v2_run'] < 100)].iterrows():
    print(row['dataset'], row['v1_project'], row['v1_run'], row['v2_project'], row['v2_run'])
```

if you have other questions or suggestions to add here, ask Matt!
