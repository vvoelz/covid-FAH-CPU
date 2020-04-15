import os, sys, glob
import numpy as np

import xvg_tools
import argparse, textwrap

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\

    Computes an MBAR estimate of the free energy of *adding* the harmonic restraint, $\Delta G_rest

    EXAMPLE
    $ python harmonic_rest_correction.py ~/server2/projects/p14602 ~/server2/data/SVR2616698070 14399 0 testout

      ''' ))

parser.add_argument('setupdir', type=str, help='The setup dir in which RUN0etup')
parser.add_argument('datadir', type=str, help='The direcory where projdata is saved')
parser.add_argument('projnum', type=int, help='The project number')
parser.add_argument('run', type=int, help='The run number')
parser.add_argument('outdir', type=str, help='The directory where outputfiles will be saved.')
parser.add_argument('--verbose', dest='verbose', action='store_true',
                    help='If specified, print output verbosely')
args = parser.parse_args()

# process the input arguments 
print('args.setupdir', args.setupdir)
index_file = os.path.join(args.setupdir, 'RUN%d/index.ndx'%args.projnum)
print('\tindex_file =', index_file)

print('args.datadir', args.datadir)
print('args.projnum', args.projnum)
print('args.run', args.run)
rundatadir = os.path.join(args.datadir, 'PROJ%d/RUN%s'%(args.projnum,args.run))
print('\trundatadir =', rundatadir)

print('args.outdir', args.outdir)
if not os.path.exists(args.outdir):
    print('Directory', args.outdir, 'does not exist.  Creating...', end='')
    os.mkdir(args.outdir)
    print('Done.')

if len(sys.argv) < 4:
    print(usage)
    sys.exit(1)


print()
print('PROJECT', args.projnum)

if (1):

    clonedirs = glob.glob( os.path.join(rundatadir, 'CLONE*') )
    clones = [ int(os.path.basename(clonedir).replace('CLONE','')) for clonedir in clonedirs]

    # sort those clones!
    Isort = np.argsort(np.array(clones))
    clonedirs = [ clonedirs[j] for j in Isort]
    clones = [ clones[j] for j in Isort ]
    print('clones', clones) 

    # for each clone, grab the data in the dhdl.xvg and pullx.xvg in each gen
    ## for j in range(len(clones)):
    
    # TESTING
    for j in [0]:
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

        # remove the last results directory, which is empty!
        resultdirs.pop() 

        print('resultdirs')
        for resultdir in resultdirs:
            print('\t',resultdir)

        ### Scrape and collect all the dhdl energies!
        for resultdir in resultdirs[0:1]:
            dhdl_xvgfile =  os.path.join(resultdir, 'dhdl.xvg')
            time_in_ps, states, energies = xvg_tools.get_dhdl(dhdl_xvgfile)
            if (args.verbose):
                print(resultdir)
                print('\ttime_in_ps', time_in_ps)
                print('\tstates', states)
                print('\tenergies.shape', energies.shape)
                print('\tenergies', energies)
            

        for resultdir in resultdirs[1:]: 
            dhdl_xvgfile =  os.path.join(resultdir, 'dhdl.xvg')
            more_time_in_ps, more_states, more_energies = xvg_tools.get_dhdl(dhdl_xvgfile)
            time_in_ps = np.concatenate( (time_in_ps, more_time_in_ps[1:]), axis=0)
            states     = np.concatenate( (states, more_states[1:]), axis=0 )
            energies   = np.concatenate( (energies, more_energies[1:,:]), axis=0 )

        states_outfile = os.path.join(args.outdir, 'r%d-c%d.states.npy'%(run,clone))
        np.save(states_outfile, states)
        print('Wrote:', states_outfile)

        energies_outfile = os.path.join(args.outdir, 'r%d-c%d.energies.npy'%(run,clone))
        np.save(energies_outfile, energies)
        print('Wrote:', energies_outfile)

        ### Scrape and collect all the pullx distances!
        for resultdir in resultdirs[0:1]:
            pullx_xvgfile =  os.path.join(resultdir, 'pullx.xvg')
            times, distances = xvg_tools.get_distances(pullx_xvgfile)
            print('distances', distances)

        for resultdir in resultdirs[1:]:
            pullx_xvgfile =  os.path.join(resultdir, 'pullx.xvg')
            more_times, more_distances = xvg_tools.get_distances(pullx_xvgfile)
            times = np.concatenate( (times, more_times[1:]), axis=0)
            distances = np.concatenate( (distances, more_distances[1:]), axis=0)

        distances_outfile = os.path.join(args.outdir, 'r%d-c%d.distances.npy'%(run,clone))
        np.save(distances_outfile, distances)
        print('Wrote:', distances_outfile)

"""
[ Restraint-Distance ]
      2527 5408
"""


    #    dhdl_xvgfiles = os.path.join(clonedir
    #try:
        
    



