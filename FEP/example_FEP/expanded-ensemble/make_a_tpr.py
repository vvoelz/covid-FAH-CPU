import os, sys, stat
from expanded import *

###############
## Let's see if grompp can build an expanded ensemble simulaton

os.system('gmx --version')

##################

# Gather all the setup files in this dir
grofile = 'npt.gro'
topfile = 'topol.top'
ndxfile = 'index.ndx'

# use the expanded_ensemble_mdpfile() class to create an initial mdpfile
e = expanded_ensemble_mdpfile()
mdpfile = 'ee.mdp'
e.write_to_filename(mdpfile)


#### make the first round of simulation
jobdir = 'out1'
if not os.path.exists(jobdir):
    os.mkdir(jobdir)

# grompp to make a tpr
grompp_cmd = "gmx grompp -c {grofile} -f {mdpfile} -p {topfile} -n {ndxfile} -o {jobdir}/frame0.tpr -po {jobdir}/mdout.mdp -maxwarn 1".format(grofile=grofile, mdpfile=mdpfile, topfile=topfile, ndxfile=ndxfile, jobdir=jobdir)
print('>>', grompp_cmd)
os.system(grompp_cmd)

# write a runme 
runme_file = os.path.join(jobdir, 'runme')
fout = open(runme_file, 'w')
fout.write('gmx mdrun -v -s frame0.tpr\n')
fout.close()
os.chmod(runme_file, stat.S_IRWXU)

print('Wrote:', runme_file)
