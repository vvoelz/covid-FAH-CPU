import os, sys, glob
import numpy as np

import xvg_tools



usage = """
    Usage:   scrape_data.py [project-data-dir] [outdir]  [runs]

EXAMPLE

    $ python scrape_data.py test-projectdata test-scraped-data 0 
"""

if len(sys.argv) < 3:
    print(usage)
    sys.exit(1)

# parse the input arguments
projdatadir = sys.argv[1]
outdir      = sys.argv[2]
runs        = [int(s) for s in sys.argv[3:]]

print('projdatadir =', projdatadir)
print('outdir =', outdir)
print('runs =', runs)

if not os.path.exists(outdir):
    print('Directory', outdir, 'does not exist.  Creating...', end='')
    os.mkdir(outdir)
    print('Done.')

# Find the viable clones 
rundirs = [ os.path.join(projdatadir, 'RUN%d'%run) for run in runs ]
print('rundirs', rundirs)

for i in range(len(runs)):
    rundir = rundirs[i]
    run = runs[i]

    clonedirs = glob.glob( os.path.join(rundir, 'CLONE*') )
    clones = [ int(os.path.basename(clonedir).replace('CLONE','')) for clonedir in clonedirs]
    print('clones', clones) 

    # for each clone, grab the data in the dhdl.xvg and pullx.xvg in each gen
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
        print('resultdirs', resultdirs)

        ### Scrape and collect all the dhdl energies!
        for resultdir in resultdirs[0:1]:
            dhdl_xvgfile =  os.path.join(resultdir, 'dhdl.xvg')
            time_in_ps, states, energies = xvg_tools.get_dhdl(dhdl_xvgfile)
            print('energies', energies)

        for resultdir in resultdirs[1:]:
            dhdl_xvgfile =  os.path.join(resultdir, 'dhdl.xvg')
            more_time_in_ps, more_states, more_energies = xvg_tools.get_dhdl(dhdl_xvgfile)
            time_in_ps = np.concatenate( (time_in_ps, more_time_in_ps[1:]), axis=0)
            states     = np.concatenate( (states, more_states[1:]), axis=0 )
            energies   = np.concatenate( (energies, more_energies[1:,:]), axis=0 )

        states_outfile = os.path.join(outdir, 'r%d-c%d.states.npy'%(run,clone))
        np.save(states_outfile, states)
        print('Wrote:', states_outfile)

        energies_outfile = os.path.join(outdir, 'r%d-c%d.energies.npy'%(run,clone))
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

        distances_outfile = os.path.join(outdir, 'r%d-c%d.distances.npy'%(run,clone))
        np.save(distances_outfile, distances)
        print('Wrote:', distances_outfile)


    #    dhdl_xvgfiles = os.path.join(clonedir
    #try:
        
    



