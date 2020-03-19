#! /usr/bin/env python3

import os, sys, glob
import shutil
from expanded import *

Testing = False      #True

usage = """
    Usage:   continue_a_tpr.py [prev-tprfile] [prev-jobdir] [this-jobdir] [topfile] [ndxfile] [ps to extend] [output tprfile] <ligonly>

    Add extra argument "ligonly" if this is a ligand-only simulation

EXAMPLE

    $ python continue_a_tpr.py testout/frame0.tpr testout nextout RUN0/topol.top  RUN0/index.ndx 1000  nextout/frame1.tpr
or
    $ python continue_a_tpr.py testout/frame0.tpr testout nextout RUN0/topol.top  RUN0/index.ndx 1000  nextout/frame1.tpr ligonly

"""

if len(sys.argv) < 8:
    print(usage)
    sys.exit(1)

print('sys.argv', sys.argv)

# parse the input arguments
prev_gen_tpr = sys.argv[1]

prev_jobdir = sys.argv[2]
prev_gen_trr = glob.glob( os.path.join(prev_jobdir, '*.trr') )[0]
#prev_gen_dhdl = os.path.join(prev_jobdir, 'dhdl.xvg')
prev_gen_mdlog = os.path.join(prev_jobdir, 'md.log')
prev_gen_ndx = os.path.join(prev_jobdir, 'index.ndx')

this_jobdir = sys.argv[3]

topfile      = sys.argv[4]
ndxfile      = sys.argv[5]
extend_in_ps = int(sys.argv[6])   # currently no being used
this_gen_tpr = sys.argv[7]

ligand_only = False
if len(sys.argv) > 8:
    if sys.argv[8].count('only') > 0:
        ligand_only = True



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

### chunk_lines IN THE NEW SHIRTS PROTOCOL should look like this:
"""             MC-lambda information
  Wang-Landau incrementor is:           8
  N   FEPL  CoulL   VdwL    Count   G(in kT)  dG(in kT)
  1  0.000  0.000  0.000       10    0.00000    4.00000   
  2  0.000  0.050  0.000       12    4.00000    0.00000   
  3  0.000  0.100  0.000       12    4.00000   16.00000   
  4  0.000  0.150  0.000       10   20.00000    8.00000   
  5  0.000  0.200  0.000        9   28.00000    2.00000   
  6  0.000  0.250  0.000       10   30.00000    2.00000   
  7  0.000  0.300  0.000       11   32.00000   10.00000   
  8  0.000  0.350  0.000       11   42.00000    0.00000   
  9  0.000  0.400  0.000       11   42.00000    0.00000   
 10  0.000  0.450  0.000       11   42.00000    2.00000   
 11  0.000  0.500  0.000       12   44.00000   14.00000   
 12  0.000  0.550  0.000        9   58.00000   -6.00000   
 13  0.000  0.600  0.000       11   52.00000   10.00000   
 14  0.000  0.650  0.000       11   62.00000    0.00000   
 15  0.000  0.700  0.000       11   62.00000    0.00000   
 16  0.000  0.750  0.000       11   62.00000   16.00000   
 17  0.000  0.800  0.000        9   78.00000    0.00000   
 18  0.000  0.850  0.000        9   78.00000    8.00000   
 19  0.000  0.900  0.000        8   86.00000    2.00000   
 20  0.000  0.950  0.000        9   88.00000    2.00000   
 21  0.000  1.000  0.000       10   90.00000    0.00000   
 22  0.000  1.000  0.100       10   90.00000   -2.00000   
 23  0.000  1.000  0.200        9   88.00000    8.00000   
 24  0.000  1.000  0.300        8   96.00000  -10.00000   
 25  0.000  1.000  0.400        8   86.00000    8.00000   
 26  0.000  1.000  0.450        7   94.00000  -10.00000   
 27  0.000  1.000  0.500        7   84.00000   10.00000   
 28  0.000  1.000  0.550        7   94.00000    0.00000   
 29  0.000  1.000  0.600        7   94.00000    6.00000   
 30  0.000  1.000  0.630        5  100.00000    0.00000   
 31  0.000  1.000  0.660        5  100.00000  -16.00000   
 32  0.000  1.000  0.690        7   84.00000    6.00000   
 33  0.000  1.000  0.720        5   90.00000    2.00000   
 34  0.000  1.000  0.750        6   92.00000   -8.00000   
 35  0.000  1.000  0.780        7   84.00000   -2.00000   
 36  0.000  1.000  0.810        6   82.00000    8.00000   
 37  0.000  1.000  0.840        5   90.00000   -6.00000   
 38  0.000  1.000  0.880        7   84.00000    8.00000 <<
 39  0.000  1.000  0.920        6   92.00000   -2.00000   
 40  0.000  1.000  1.000        5   90.00000    0.00000   

"""

# Discard the first three lines and pull out the current energy incrementor...
chunk_lines.pop(0)    #          MC-lambda information
wl_increment_in_kT = float(chunk_lines.pop(0).split()[-1])     #   Wang-Landau incrementor is:           3
chunk_lines.pop(0)    #  N   FEPL    Count   G(in kT)  dG(in kT)
print('current Wang Landau increment_in_kT', wl_increment_in_kT)
# ... and values from the 4th column to get the latest Wang Landau weights


nlambdas = len(chunk_lines)
## row index 5 has the WL weights                            vvv !!   
wang_landau_weights = np.array([ float(chunk_lines[i].split()[5]) for i in range(nlambdas) ])
print('current wang_landau_weights', wang_landau_weights)


# Write an mdp file with the latest weights !!!
e = expanded_ensemble_mdpfile(ligand_only=ligand_only, init_lambda_weights=wang_landau_weights,
                              wl_increment_in_kT=wl_increment_in_kT)
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

