import os, sys, glob
import numpy as np

import xvg_tools


import argparse, textwrap

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\

    Shows progress bar of how much sampling the simulations have reached

    $ python progress_bar.py  /home/server/server2/data/SVR2616698070  14601 

      ''' ))

parser.add_argument('datadir', type=str, help='the direcory where projdata is saved')
parser.add_argument('projnum', type=int, help='the project number')
#parser.add_argument('--agg', dest='agg', action='store_true',
#                    help='Specify the structure comes from an AGG simulation.  The resnum in conf.gro is off by +1!')

args = parser.parse_args()
print('args.datadir', args.datadir)
print('args.projnum', args.projnum)


# find the number of runs from the project.xml   file
project_xml_file = '/home/server/server2/projects/p%d/project.xml'%args.projnum
fin = open(project_xml_file, 'r')
lines  = fin.readlines()
for line in lines:
    if line.count('runs'):
        fields = line.split('"')
        nruns = int(fields[1])

# Find the viable clones 
rundirs = [ os.path.join(args.datadir, 'PROJ%d/RUN%d'%(args.projnum,run)) for run in range(0,nruns) ]

print()
print('PROJECT', args.projnum)

for i in range(nruns):

    rundir = rundirs[i]
    clonedirs = glob.glob( os.path.join(rundir, 'CLONE*') )
    clones = [ int(os.path.basename(clonedir).replace('CLONE','')) for clonedir in clonedirs]
    run_label = 'RUN%d'%i
    #run_header = '%-8s | .....................................50|ns..................................100|ns.............................................150|ns'%run_label
    run_header = '%-8s | .....................................50|ns..................................100|ns'%run_label
    print(run_header)

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

        print('CLONE %04d'%clone, '*'*len(resultdirs))

    #    dhdl_xvgfiles = os.path.join(clonedir
    #try:
        
    



