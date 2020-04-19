#! /usr/bin/env python3

import os, sys, glob, stat, shutil
import subprocess
import errno

import numpy as np

from expanded_v3 import *

###############################################################
### First thing's first: do we have a GROMACS install or what?
global GMX_BIN
try:
    GMX_BIN = os.environ['GMXBIN']
except:
    GMX_BIN = '/usr/local/gromacs/bin'
if not os.path.exists(GMX_BIN):
    print('Cannot find the GROMACS installation!  Please set environment variable GMX_BIN')
    sys.exit(1)
###############################################################


# Subroutines and Functions

def most_converged_mdlogfile(project_run_datadir):
    """Finds and returns the md.log file with the most converged expanded-ensemble sampling 
    By "most converged", we mean the longest number of gens. If there is a tie, we just pick the first one
    encountered.

    INPUT
    project_run_datadir      - pathname of ~/server2/data/SVRxxxxxxx/PROJxxxx/RUNxxx

    OUTPUT
    mdlogfile                - pathname of the md.log file
    """

    ###################################################
    # Let's find the CLONE directory with largest gen!

    clonedirs = glob.glob( os.path.join(project_run_datadir, 'CLONE*') )
    clones = [ int(os.path.basename(clonedir).replace('CLONE','')) for clonedir in clonedirs]

    # sort those clones!
    Isort = np.argsort(np.array(clones))
    clonedirs = [ clonedirs[j] for j in Isort]
    clones = [ clones[j] for j in Isort ]

    # for each clone, count the number of results dirs
    maxgen_per_clone = []
    for j in range(len(clones)):
        clone = clones[j]
        clonedir = clonedirs[j]

        ### Find all the result0, result1, etc. dirs !
        resultdirs = []
        resultdirs1 = glob.glob( os.path.join(clonedir, 'results?') )
        resultdirs1.sort()
        resultdirs += resultdirs1

        resultdirs2 = glob.glob( os.path.join(clonedir, 'results??') )
        resultdirs2.sort()
        resultdirs += resultdirs2

        resultdirs3 = glob.glob( os.path.join(clonedir, 'results???') )
        resultdirs3.sort()
        resultdirs += resultdirs3

        # the very last resultdir is empty, because of how the continue_a_tpr_vX.py scripts work
        # so lets ignore it 
        resultdirs.pop()
        print('resultdirs', resultdirs)

        maxgen_per_clone.append( len(resultdirs)-1 ) 

    # which clone is the longest?
    maxgen_per_clone = np.array( maxgen_per_clone )
    maxgen_clone_index  = np.argmax( maxgen_per_clone )
    
    most_converged_resultdir = os.path.join(clonedirs[maxgen_clone_index], 'results%d'%maxgen_per_clone[maxgen_clone_index])
    print('most_converged_resultdir', most_converged_resultdir) 

    mdlogfile = os.path.join(most_converged_resultdir, 'md.log')
    if not os.path.exists(mdlogfile):
        print("Can't find the md.log file!", mdlogfile)
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), mdlogfile)

    return mdlogfile



def extract_wl_info_from_mdlog(mdlogfile):
    """Extracts the wl_increment_in_kT and wl_weights from a md.log file
    
    INPUT
    mdlogfile          - pathname of the input md.log file
    
    OUTPUT
    wl_increment_in_kT - (float) value of the last wl_increment in the md.log file
    wl_weights         - Wang-Landau weights as np.array()
    wl_state_index     - the last wl_lambda state index
    """

    fin = open(mdlogfile, 'r')
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
    column_headers = (chunk_lines[0].strip()).split()    #  N   FEPL    Count   G(in kT)  dG(in kT)
    print('column_headers', column_headers)
    WL_weight_column_index = column_headers.index('G(in')
    chunk_lines.pop(0)
    print('current Wang Landau increment_in_kT', wl_increment_in_kT)
    # ... and values from the 4th column to get the latest Wang Landau weights
    
    # Parse the WL weights
    nlambdas = len(chunk_lines)
    ## Using the MRS protocol, row index 5 has the WL weights    VVV   !!   
    ## if there are FEP rest lambdas, it's row 6, e.g.
    ## USE THE column index found above 
    wang_landau_weights = np.array([ float(chunk_lines[i].split()[WL_weight_column_index]) for i in range(nlambdas) ])
    print('current wang_landau_weights', wang_landau_weights)

    # Get the latest lambda state index
    wl_state_index = 0
    for line in chunk_lines:
        if line.count( ' <<' ) > 0:
            wl_state_index = int(line.split()[0]) - 1  # Argh!!!  the numbering started at 1!  Fixed now

    return wl_increment_in_kT, wang_landau_weights, wl_state_index 

#########################

# Main

if __name__ == '__main__':

    import argparse, textwrap

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\

    init_wl_weights.py

    This script will take an input *.mdp file and replace (or add) the initial Wang-Landau weights with the
    longest-running / most converged weights found across the clones in the specified ...data/SVRxxxx/PROJxxx/RUNxxx directory

    If input *.mdp is the same as the output *.mdp, a backup *.mdp -> *.mdp.backup will be written 
    
    EXAMPLE
    $ python init_wl_weights.py RUN0/prod.mdp /home/server/server2/data/SVR2616698070/PROJ14399/RUN0 RUN0/prod.mdp

      ''' ))

    parser.add_argument('input_mdpfile', type=str, help='the input *.mdp file')
    parser.add_argument('project_run_datadir', type=str, help='The FAH data directory e.g. ~/server2/data/SVRxxxxxxx/PROJxxxx/RUNxxx')
    parser.add_argument('output_mdpfile', type=str, help='the output *.mdp file')

    args = parser.parse_args()
    print('args.input_mdpfile', args.input_mdpfile)
    print('args.project_run_datadir', args.project_run_datadir)
    print('args.output_mdpfile', args.output_mdpfile)

    # Find the most converged mdlogfile
    mdlogfile = most_converged_mdlogfile(args.project_run_datadir)

    # get the WL info from it
    wl_increment_in_kT, wang_landau_weights, wl_state_index = extract_wl_info_from_mdlog(mdlogfile)

    # read in the old mdpfile
    e = expanded_ensemble_mdpfile()
    e.read_parms_from_mdpfile(args.input_mdpfile)
    #print('##### input_mdpfile', args.input_mdpfile, '#####')
    #e.report()

    # store the new WL settings
    e.init_lambda_weights  = wang_landau_weights
    e.init_lambda_weights_string =  ' '.join(['%2.5f'%e.init_lambda_weights[i] for i in range(e.nlambdas)])
    e.init_lambda_state    = wl_state_index
    e.wl_increment_in_kT   = wl_increment_in_kT

    # write the new mdpfile
    if os.path.exists(args.output_mdpfile):
        print('Output *.mdp file', args.output_mdpfile, 'already exists!')
        backup_mdpfile = args.output_mdpfile+'.backup'
        os.system('mv {output_mdpfile} {backup_mdpfile}' .format(output_mdpfile=args.output_mdpfile, backup_mdpfile=backup_mdpfile) )
        print('...moved to:', backup_mdpfile)

    e.write_to_filename(args.output_mdpfile)
    print()
    print('*** Wrote:', args.output_mdpfile)
    
