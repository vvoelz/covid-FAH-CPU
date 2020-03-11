#!/usr/bin/env

import os, shutil

for i in range(1,8):
    os.mkdir(f'HIT{i}')
    for file in [f'hit{i}.sdf',f'hit{i}.mol2',f'hit{i}.pdb','receptor.mol2','receptor.pdb']:
        shutil.copy2(file,f'HIT{i}')
