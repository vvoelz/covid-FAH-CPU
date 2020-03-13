#! /usr/bin/env python3

import os, sys
import shutil
from expanded import *

Testing = False      #True

usage = """
    Usage:   continue_a_tpr.py [prev-jobdir] [this-jobdir] [ndxfile] [ps to extend]

EXAMPLE

    $ python continue_a_tpr.py out1 out2 index.ndx 1000
"""

if len(sys.argv) < 5:
    print(usage)
    sys.exit(1)

# parse the input arguments
prev_jobdir = sys.argv[1]
prev_gen_tpr = os.path.join(prev_jobdir, 'frame0.tpr')   # <-- we assume a certain filename convention
prev_gen_trr = os.path.join(prev_jobdir, 'traj.trr')
#prev_gen_dhdl = os.path.join(prev_jobdir, 'dhdl.xvg')
prev_gen_mdlog = os.path.join(prev_jobdir, 'md.log')
prev_gen_ndx = os.path.join(prev_jobdir, 'index.ndx')

this_jobdir = sys.argv[2]
this_gen_tpr = os.path.join(this_jobdir, 'frame1.tpr')

ndxfile      = sys.argv[3]
extend_in_ps = int(sys.argv[4])

#############################
# Read in the previous md.log file and parse out the latest WL weights!

fin = open(prev_gen_mdlog, 'r')
lines = fin.readlines()
fin.close()

# find the lines with the string "MC-lambda information"
start_indices = [i for i in range(len(lines)) if lines[i].count("MC-lambda information") > 0 ]
j = start_indices[-1] # starting index of the last chunk
chunk_lines = []
while lines[j].strip() != '':
    chunk_lines.append(lines[j])
    j += 1
print(''.join(chunk_lines))

### chunk_lines should look like this:
"""          MC-lambda information
  Wang-Landau incrementor is:           3
  N   FEPL    Count   G(in kT)  dG(in kT)
  1  0.000       24    0.00000    0.00000   
  2  0.050       24    0.00000   -3.00000   
  3  0.100       25   -3.00000    0.00000 <<
  4  0.150       25   -3.00000   -6.00000   
  5  0.200       27   -9.00000   -9.00000   
  6  0.250       30  -18.00000   -3.00000   
  7  0.300       31  -21.00000   -9.00000   
  8  0.350       34  -30.00000  -12.00000   
  9  0.400       38  -42.00000   -6.00000   
 10  0.450       40  -48.00000  -12.00000   
 11  0.500       44  -60.00000  -12.00000   
 12  0.550       48  -72.00000  -12.00000   
 13  0.600       52  -84.00000   -9.00000   
 14  0.650       55  -93.00000   -9.00000   
 15  0.700       58 -102.00000  -15.00000   
 16  0.750       63 -117.00000   -9.00000   
 17  0.800       66 -126.00000  -18.00000   
 18  0.850       72 -144.00000  -12.00000   
 19  0.900       76 -156.00000  -15.00000   
 20  0.950       81 -171.00000  -18.00000   
 21  1.000       87 -189.00000    0.00000   
"""
# Discard the first three lines and pull out the current energy incrementor...
chunk_lines.pop(0)    #          MC-lambda information
wl_increment_in_kT = float(chunk_lines.pop(0).split()[-1])     #   Wang-Landau incrementor is:           3
chunk_lines.pop(0)    #  N   FEPL    Count   G(in kT)  dG(in kT)
print('current Wang Landau increment_in_kT', wl_increment_in_kT)
# ... and values from the 4th column to get the latest Wang Landau weights


nlambdas = len(chunk_lines)
wang_landau_weights = np.array([ float(chunk_lines[i].split()[3]) for i in range(nlambdas) ])
print('current wang_landau_weights', wang_landau_weights)

# write an mdp file with the latest weights
e = expanded_ensemble_mdpfile(init_lambda_weights=wang_landau_weights,
                              wl_increment_in_kT=wl_increment_in_kT)
if not os.path.exists(this_jobdir):
    os.mkdir(this_jobdir)
this_mdpfile = os.path.join(this_jobdir, 'ee.mdp') 
print('Writing', this_mdpfile, '...')
e.write_to_filename(this_mdpfile)
print ('...Done')


cmd = "gmx grompp -c {prev_gen_tpr} -t {prev_gen_trr} -f {this_mdpfile} -n {ndxfile} -o {this_gen_tpr} -po {this_jobdir}/mdout.mdp -maxwarn 1".format(prev_gen_tpr=prev_gen_tpr, prev_gen_trr=prev_gen_trr, \
	this_mdpfile=this_mdpfile, ndxfile=ndxfile, this_gen_tpr=this_gen_tpr, this_jobdir=this_jobdir)
print('cmd:', cmd)
# os.system(cmd)

