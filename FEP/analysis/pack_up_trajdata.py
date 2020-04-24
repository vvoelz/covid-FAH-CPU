import os, sys, glob
import numpy as np

import xvg_tools
import traj_tools
import subprocess

import argparse, textwrap


global gmx_bin
gmx_bin = 'gmx'

### Functions ###

def run_cmd(cmd, testing=False, verbose=True):
    """Print and do an os.system() on the input command"""

    if verbose:
        print('>>', cmd)
    if not testing:
        subprocess.check_output(cmd, shell=True)

def check_and_mkdir(dirpath, testing=False, verbose=True):
    """Check to see of the directory exits.
    If not, use os.mkdir() to make a new directory, with fancy verbosity."""

    print('os.path.exists(dirpath)', os.path.exists(dirpath))

    if not os.path.exists(dirpath):
        print('Directory', dirpath, 'does not exist...', end='')
        os.mkdir(dirpath)
        print('Created directory', dirpath, '\n')


def run_cmd_interactive(cmd, input_string, testing=False, verbose=True):
    """Run the command and pass it an interactive input string.

    NOTES
    from https://stackoverflow.com/questions/163542/how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument

    EXAMPLE
    p = subprocess.Popen(['grep', 'f'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    grep_stdout = p.communicate(input=b'one\ntwo\nthree\nfour\nfive\nsix\n')[0]
    print(grep_stdout.decode())
    """

    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    grep_stdout = p.communicate(input=input_string.encode('utf-8'))[0]
    if verbose:
        print(grep_stdout.decode())





### Main ####

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\

    Packs up all the relevant trajectory data for Protein + LIG system (stripping HOH, NA, CL)
    Files written to the outdir:

    outdir
    | 
    |---RUN0
        |---c0.distances.npy
        |---c0.energies.npy 
        |---c0.states.npy   
        |---c0.xtc
        |
        |---c1.distances.npy
        |---c1.energies.npy
        |---c1.states.npy
        |---c1.xtc
        |
        | ... 
        |
        |---c9.distances.npy
        |---c9.energies.npy
        |---c9.states.npy
        |---c9.xtc
        |
        |---Protein_LIG.tpr
        |---Protein_LIG.ndx
        |---Protein_LIG.gro
         
    |---RUN1
        |---c0.distances.npy
        ...
    ...

    EXAMPLE
    $ python pack_up_trajdata.py /home/server/server2/projects/p14399 ~/server2/data/SVR2616698070 14399 -r 0,1,3-4 p14399-trajdata

      ''' ))

parser.add_argument('setupdir', type=str, help='The setup directory for this project, e.g. ~/server2/projects/p14399')
parser.add_argument('datadir', type=str, help='The direcory where projdata is saved')
parser.add_argument('projnum', type=int, help='The project number')
#parser.add_argument('--agg', dest='agg', action='store_true',
#                    help='Specify the structure comes from an AGG simulation.  The resnum in conf.gro is off by +1!')
parser.add_argument('-r', action='append',
                    dest='specific_runs',
                    default=[],
                    help='A string (with no spaces! e.g. 0,1,3-4) denoting a specific set of runs')
parser.add_argument('outdir', type=str, help='The direcory where all the projdata is saved')
parser.add_argument('--verbose', dest='verbose', action='store_true',
                    help='If specified, print output verbosely')


args = parser.parse_args()
if args.verbose:
    print('### Input arguments parsed ###')
    print('args.setupdir', args.setupdir)
    print('args.datadir', args.datadir)
    print('args.projnum', args.projnum)
    print('args.specific_runs', args.specific_runs)
    print('args.outdir', args.outdir)
    print('\n')

# If the output directory doesn't exist, make it
check_and_mkdir(args.outdir)


## If the user has chosen specific runs, parse this string
myruns = None
if len(args.specific_runs) > 0:
    use_specific_runs = True
    myruns = []
    fields = args.specific_runs[0].split(',')
    for field in fields:
        if field.count('-') == 0:
            myruns.append(int(field))
        else:
            ends = [int(s) for s in field.split('-')]
            myruns += list(range(ends[0],ends[1]+1))
    print('User has selected to pack up runs:', myruns, '\n')

else:
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
    myruns = range(0,nruns)
    print('Packing up all runs:', myruns, '\n') 

rundirs = [ os.path.join(args.datadir, 'PROJ%d/RUN%d'%(args.projnum,run)) for run in myruns ]

for i in range(len(myruns)):

    run = myruns[i]
    rundatadir = rundirs[i]

    this_run_outdir = os.path.join(args.outdir, 'RUN%d'%run)
    check_and_mkdir(this_run_outdir)

    # write the Protein_LIG.tpr to the output directory  
    this_setup_rundir = os.path.join(args.setupdir, 'RUN%d'%run)
    traj_tools.build_Protein_LIG_files(this_setup_rundir, this_run_outdir)


    clonedirs = glob.glob( os.path.join(rundatadir, 'CLONE*') )
    clones = [ int(os.path.basename(clonedir).replace('CLONE','')) for clonedir in clonedirs]

    # sort those clones!
    Isort = np.argsort(np.array(clones))
    clonedirs = [ clonedirs[j] for j in Isort]
    clones = [ clones[j] for j in Isort ]
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

        
        # remove the last results directory, which is empty!
        if len(resultdirs) > 0:
            resultdirs.pop() 

        # If --lastgenonly, only look at the last resultsdir
        #if args.lastgenonly:
        #    resultdirs = [resultdirs[-1]]

        print('resultdirs')
        for resultdir in resultdirs:
            print('\t',resultdir)

        ###### XTC ########

        ### concatenate all the xtc files into one!
        out_xtcfile = os.path.join(this_run_outdir, 'c%d.xtc'%clone)

        # make a list of all the xtc_files:
        result_xtcfiles = [os.path.join(resultdir,'traj_comp.xtc') for resultdir in resultdirs if os.path.exists(os.path.join(resultdir,'traj_comp.xtc'))]
        if len(result_xtcfiles) > 0: 
            cmd = '{gmx} trjcat -o {out_xtcfile} -f '
            for result_xtcfile in result_xtcfiles:
                cmd += (result_xtcfile + ' ')
            cmd += '-cat'
            run_cmd(cmd.format(gmx=gmx_bin, out_xtcfile=out_xtcfile))
            # fix PBC artifacts with cluster
            out_cluster_xtcfile = os.path.join(this_run_outdir, 'c%d.cluster.xtc'%clone)
            tprfile = os.path.join(this_run_outdir, 'Protein_LIG.tpr')
            cmd = '{gmx} trjconv -s {tprfile} -o {out_cluster_xtcfile} -f {out_xtcfile} -pbc cluster'.format(
                    gmx=gmx_bin, tprfile=tprfile, out_cluster_xtcfile=out_cluster_xtcfile, out_xtcfile=out_xtcfile)
            run_cmd_interactive(cmd, "0\n0\n")

            # while we're at it, let's fix the PBC artifact with Protein_LIG.gro as well
            grofile = os.path.join(this_run_outdir, 'Protein_LIG.gro') 
            out_grofile = os.path.join(this_run_outdir, 'Protein_LIG.cluster.gro')        
            cmd = '{gmx} trjconv -s {tprfile} -o {out_grofile} -f {grofile} -pbc cluster'.format(
                gmx=gmx_bin, tprfile=tprfile, out_grofile=out_grofile, grofile=grofile)
            run_cmd_interactive(cmd, "0\n0\n")

        ###### ENERGIES ########

        result_dhdlfiles = [os.path.join(resultdir,'dhdl.xvg') for resultdir in resultdirs if os.path.exists(os.path.join(resultdir,'dhdl.xvg'))]
        if len(result_dhdlfiles) > 0:

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

            states_outfile = os.path.join(this_run_outdir, 'c%d.states.npy'%clone)
            np.save(states_outfile, states)
            print('Wrote:', states_outfile)

            energies_outfile = os.path.join(this_run_outdir, 'c%d.energies.npy'%clone)
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

            distances_outfile = os.path.join(this_run_outdir, 'c%d.distances.npy'%clone)
            np.save(distances_outfile, distances)
            print('Wrote:', distances_outfile)
    
### Cleanup

os.system('rm -f ./#mdout*')



