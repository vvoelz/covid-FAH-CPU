#! /usr/bin/env python3

import os, sys, glob
import shutil

from init_wl_weights import *
from expanded_v2 import *

Testing = False      #True

usage = """
    Usage:   continue_a_tpr_v2.py [prev-tprfile] [prev-jobdir] [this-jobdir] [topfile] [ndxfile] [mdpfile] [output tprfile] <ligonly>

    Add extra argument "ligonly" if this is a ligand-only simulation

    VERSION NOTES
    04-13-2020   v2: This version of the script, and the helper object expanded_v2.py are designed to FIX problems with the umbrella restraint.
                     NOTE the reliance on init_wl_weights.py methods to extract from WL info from md.log files 

    EXAMPLE
        $ python continue_a_tpr_v2.py runclone0/frame0.tpr runclone0/results0 nextout RUN0/topol.top  RUN0/index.ndx RUN0/prod.mdp  nextout/frame1.tpr
    or
        $ python continue_a_tpr_v2.py runclone0/frame0.tpr runclone0/results0 nextout RUN0/topol.top  RUN0/index.ndx RUN0/prod.mdp  nextout/frame1.tpr ligonly

"""

if len(sys.argv) < 8:
    print(usage)
    sys.exit(1)

print('sys.argv', sys.argv)

# parse the input arguments
prev_gen_tpr = sys.argv[1]

prev_jobdir = sys.argv[2]
prev_gen_trr = glob.glob( os.path.join(prev_jobdir, '*.trr') )[0]
prev_gen_mdlog = os.path.join(prev_jobdir, 'md.log')
prev_gen_ndx = os.path.join(prev_jobdir, 'index.ndx')

this_jobdir = sys.argv[3]

topfile      = sys.argv[4]
ndxfile      = sys.argv[5]
mdpfile      = sys.argv[6]
this_gen_tpr = sys.argv[7]

ligand_only = False
if len(sys.argv) > 8:
    if sys.argv[8].count('only') > 0:
        ligand_only = True


#############################
# Read in the previous md.log file and parse out the latest WL weights!
wl_increment_in_kT, wang_landau_weights, init_lambda_state = extract_wl_info_from_mdlog(prev_gen_mdlog)

### new in v2: get the pull_coord1_k  and pull_coord1_init !
e_from_setup = expanded_ensemble_mdpfile( ligand_only=ligand_only )
e_from_setup.read_parms_from_mdpfile(mdpfile)
use_this_pull_coord1_k    = e_from_setup.pull_coord1_k
use_this_pull_coord1_init = e_from_setup.pull_coord1_init
print('Will use pull_coord1_k =', use_this_pull_coord1_k)
print('Will use pull_coord1_init =', use_this_pull_coord1_init)

# Write an mdp file with the latest weights !!!
e = expanded_ensemble_mdpfile( ligand_only=ligand_only,
                               init_lambda_weights=wang_landau_weights,
                               init_lambda_state=init_lambda_state,
                               wl_increment_in_kT=wl_increment_in_kT,
                               pull_coord1_k = use_this_pull_coord1_k,
                               pull_coord1_init = use_this_pull_coord1_init) 

if not os.path.exists(this_jobdir):
    os.mkdir(this_jobdir)
this_mdpfile = os.path.join(this_jobdir, 'ee.mdp') 
print('Writing', this_mdpfile, '...')
e.write_to_filename(this_mdpfile)
print ('...Done')

GMX_BIN = '/usr/local/gromacs/bin/gmx'
cmd = "{GMX_BIN} grompp -c {prev_gen_tpr} -t {prev_gen_trr} -f {this_mdpfile} -p {topfile} -n {ndxfile} -o {this_gen_tpr} -po {this_jobdir}/mdout.mdp -maxwarn 1".format(GMX_BIN=GMX_BIN, prev_gen_tpr=prev_gen_tpr, prev_gen_trr=prev_gen_trr, \
	this_mdpfile=this_mdpfile, topfile=topfile, ndxfile=ndxfile, this_gen_tpr=this_gen_tpr, this_jobdir=this_jobdir)
print('cmd:', cmd)
os.system(cmd)

