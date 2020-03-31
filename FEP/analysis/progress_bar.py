import os, sys, glob
import numpy as np

import xvg_tools
import argparse, textwrap

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\

    Shows progress bar of how much sampling the simulations have reached

    EXAMPLE
    $ python progress_bar.py ~/server2/data/SVR2616698070 14602

      ''' ))

parser.add_argument('datadir', type=str, help='the direcory where projdata is saved')
parser.add_argument('projnum', type=int, help='the project number')
#parser.add_argument('--agg', dest='agg', action='store_true',
#                    help='Specify the structure comes from an AGG simulation.  The resnum in conf.gro is off by +1!')

args = parser.parse_args()
print('args.datadir', args.datadir)
print('args.projnum', args.projnum)




def project_length_in_ns(projnum):
    """Returns a float with the project length in nanoseconds (ns)"""

    w = {}

    # 7 fragment hits
    w[14346] = 1.0    # RL
    w[14348] = 10.0   # L

    #  100 ligands
    w[14337] = 1.0    # RL
    w[14339] = 10.0   # L

    # 72 series
    for i in range(14600, 14613):
        w[i] = 1.0    # RL
    for i in range(14349, 14362):
        w[i] = 10.0   # L

    # Moonshot 03-23
    w[14363] = 1.0    # RL
    w[14364] = 10.0   # L

    # Moonshot 03-26
    for i in range(14365, 14369):
        w[i] = 1.0    # RL
    for i in range(14369, 14373):
        w[i] = 10.0   # L

    # MLTN
    w[14373] = 1.0    # RL
    w[14374] = 10.0   # L

    return w[projnum]



# find the number of runs from the project.xml   file
project_xml_file = '/home/server/server2/projects/p%d/project.xml'%args.projnum
if not os.path.exists(project_xml_file):
    # ns338286 projects do not have a leading "p" .... try that
    project_xml_file = '/home/server/server2/projects/%d/project.xml'%args.projnum
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

        print('CLONE %04d'%clone, '*'*((len(resultdirs)-1)*int(project_length_in_ns(args.projnum))) )  # added -1 to account for empty results dir at end

    #    dhdl_xvgfiles = os.path.join(clonedir
    #try:
        
    



