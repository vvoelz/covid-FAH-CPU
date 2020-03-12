import os, sys
import subprocess

###############
## Let's see if grompp can build an expanded ensemble simulaton

os.system('gmx --version')

home  = os.path.curdir
jobdir = 'out1'
if not os.path.exists(jobdir):
    os.mkdir(jobdir)

grofile = 'npt.gro'
mdpfile = 'ee-umbrella.mdp'
topfile = 'topol.top'
ndxfile = 'index.ndx'


grompp_cmd = "gmx grompp -c {grofile} -f {mdpfile} -p {topfile} -n {ndxfile} -o {jobdir}/frame0.tpr -po {jobdir}/mdout.mdp -maxwarn 1".format(grofile=grofile, mdpfile=mdpfile, topfile=topfile, ndxfile=ndxfile, jobdir=jobdir)
print('grompp_cmd:')
os.system(grompp_cmd)
