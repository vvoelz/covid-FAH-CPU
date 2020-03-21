#! /usr/bin/env python

import os, sys
import subprocess, time

TESTING = False   # True


def run_cmd(cmd, testing=True, output_file=None, error_file=None):
    """Run the command. (or if testing == True, print the command.
    
    INPUT
    cmd             - a string containing the commmand

    OPTIONS
    output_file     -  if filename specified, write std output to that file
    error_file      -  if filename specified, write std error to that file
    
    RETURNS
    returncode      - 0 is successful, 1 if errors
    output          - standard output text
    errors          - standard error text
    time_in_seconds - the elapsed time in seconds
    
    """

    print('>>', cmd)
    print('>> Writing stdout to: %s ; stderr to: %s ...'%(output_file, error_file))

    if output_file == None:
        output_file = '/dev/null'
    if error_file == None:
        error_file = '/dev/null'

    if not testing:
        #os.system(cmd)
        args = cmd.split()
        #c = subprocess.run(args, capture_output=True, text=True)
        #proc = subprocess.Popen([sys.executable, cmd], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        #stdout, stderr = proc.communicate()

        start_time = time.time()

        # proc.returncode, stdout, stderr
        with open(output_file,'w+') as fout:
            with open(error_file,'w+') as ferr:
                returncode = subprocess.call(args, stdout=fout, stderr=ferr)
                # reset file to read from it
                fout.seek(0)
                # save output (if any) in variable
                output=fout.read()

                # reset file to read from it
                ferr.seek(0)
                # save errors (if any) in variable
                errors = ferr.read()

    print('>> ...Done.')
    time_in_seconds = (time.time() - start_time)
    print(">> Took --- %s seconds ---" % (time.time() - start_time))

    return returncode, output, errors, time_in_seconds


############################################

outready_dir = 'out-ready'
outcreate_dir = 'out-create'
if not os.path.exists(outready_dir):
    os.mkdir('out-ready')
if not os.path.exists(outcreate_dir):
    os.mkdir('out-create')

ready_path = '/home/server/git/covid-FAH-CPU/FEP/scripts/make_fep_ready_argparse.py'
#ready_path = '/home/server/git/covid-FAH-CPU/FEP/scripts/make_fep_ready.py'

create_path = '/home/server/git/covid-FAH-CPU/FEP/scripts/create_ee_mdp.py'

nruns = 95 

if (0):
    runs = range(10, 20)
    babysteps = False
else:
    ## Refine the ones that didn't work!
    runs = [15]
    babysteps = True

for r in runs:

    # Let's make the p14600 project FEP-ready
    if babysteps:
        cmd = 'python {ready_path} /data/72/projects/p14600/RUN{r} ./RUN{r} --babysteps'.format(ready_path=ready_path, r=r)
    else:
        cmd = 'python {ready_path} /data/72/projects/p14600/RUN{r} ./RUN{r}'.format(ready_path=ready_path, r=r)

    returncode, output, errors, time_in_seconds = run_cmd(cmd, testing=TESTING, output_file=os.path.join(outready_dir,'out.%d'%r), error_file=os.path.join(outready_dir,'err.%d'%r))
    print('RUN%d took %3.3f seconds.'%(r,time_in_seconds))

    # First, we create a protein-ligand *.mdp for the *.tpr we're about to build
    cmd = 'python {create_path} RUN{r}/npt.gro RUN{r}/topol.top RUN{r}/index.ndx RUN{r}/prod.mdp'.format(create_path=create_path, r=r)
    run_cmd(cmd, testing=TESTING, output_file=os.path.join(outcreate_dir,'out.%d'%r), error_file=os.path.join(outcreate_dir,'err.%d'%r))

    # Cleanup
    os.system('rm ./#* ; rm step* ; rm traj* ; rm ener.edr')

    print()


