import os, sys, stat
import subprocess

from expanded_v2 import *


usage = """
    Usage:   create_ee_mdp.py [grofile] [topfile] [ndxfile] [output mdpfile] ligonly

EXAMPLE

    $ python create_ee_mdp.py RUN0/npt.gro RUN0/topol.top RUN0/index.ndx RUN0/prod.mdp
or
    $ python create_ee_mdp.py RUN0/npt.gro RUN0/topol.top RUN0/index.ndx RUN0/prod.mdp ligonly
"""

if len(sys.argv) < 5:
    print(usage)
    sys.exit(1)

# parse the input arguments
prev_jobdir = sys.argv[1]
# Gather all the setup files in this dir
grofile = sys.argv[1]
topfile = sys.argv[2]
ndxfile = sys.argv[3]
mdpfile = sys.argv[4]

ligand_only = False
if len(sys.argv) > 5:
    if sys.argv[5].count('only') > 0:
        ligand_only = True



if not ligand_only:

    # To tether the ligand to the protein, we need to have prepared an index file with
    # the following atom groups, for example:
    """
    [ a1-Protein ]
    678
    [ a2-Ligand ]
    1564
    [ Restraint-Distance ]
    678 1564
    """
    # NOTE that for our tool chain to work smoothly, we should always use these DEFAULT 
    # atom group names

    ### To create the mdpfile, we need to know these default names, *and* the distance between the atoms
    GMX_BIN = '/usr/local/gromacs/bin/gmx'
    if not os.path.exists(GMX_BIN):
        GMX_BIN = 'gmx'    
    
    # write a temp "selection.dat"
    fout = open('selection.dat', 'w')
    fout.write('group "Restraint-Distance"\n')
    fout.close()
    distance_cmd = '{GMX_BIN} distance -f {grofile} -n {ndxfile} -sf selection.dat'.format(GMX_BIN=GMX_BIN, grofile=grofile, ndxfile=ndxfile)
    output = subprocess.check_output(distance_cmd, shell=True).decode("utf-8")
    output_lines = output.split('\n')
    atom_distance = None
    for line in output_lines:
        if line.count('Average') > 0:
            fields = line.split()
            print('fields', fields)
            atom_distance = float( fields[2] )  # ['Average', 'distance:', '0.387', 'nm']
    print('atom_distance is', atom_distance, 'nm')
    
    # use the expanded_ensemble_mdpfile() class to create an initial mdpfile
    e = expanded_ensemble_mdpfile(ligand_only = ligand_only,
                        pull_group1_name = 'a1-Protein',
                        pull_group2_name = 'a2-Ligand',
                        pull_coord1_init  = atom_distance)

else:
    e = expanded_ensemble_mdpfile(ligand_only = ligand_only)

e.write_to_filename(mdpfile)

