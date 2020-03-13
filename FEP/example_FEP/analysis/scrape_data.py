import os, sys, glob
import numpy as np

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

    # for each clone, grab the data in the dhdl.xvg and pullfind the 
    



