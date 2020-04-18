import os, sys, glob
import errno
import doctest
import numpy as np

from time import sleep
import subprocess


######## Functions #############

def run_cmd(cmd, testing=False):
    """Print and do an os.system() on the input command"""

    print('>>', cmd)
    if not testing:
        subprocess.check_output(cmd, shell=True)


def build_Protein_LIG_tpr(setup_rundir, out_tprfile, version='v1', gmx_bin='gmx', verbose=False):
    """Builds a Protein_LIG only tpr for the purpose of looking at trajectories
    in Chimera, PBC transformations, etc.
    
    INPUT
    setup_rundir - a setup RUN directory like /home/server/server2/projects/p14399/RUN0/
    out_tprfile  - filename to write the tpr file

    PARAMETERS
    version      - 'v1' for the original FEP mdps, 'v2' for the _v2 protocol, etc. (Default: 'v1')
    gmx_bin      - pathname to the gmx executable (Default: 'gmx')
    verbose      - if True , print verbose output

    USAGE

    >>> build_Protein_LIG_tpr('/home/server/server2/projects/p14399/RUN0')
    """

    ### Create a new index group for Protein_LIG ##

    # find the ndxfile
    ndxfile = os.path.join(setup_rundir, 'index.ndx')
    if not os.path.exists(ndxfile):
        print("Can't find index.ndx file", ndxfile)
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), ndxfile)

    # Get the indices of Protein and LIG
    index_Protein, index_LIG = get_ndx_group_index(ndxfile, ['Protein', 'LIG'])
    print('index_Protein', index_Protein, 'index_LIG', index_LIG)

    # make a new atom group that is the union of Protein nd LIG
    out_ndxfile = os.path.join(setup_rundir, 'Protein_LIG.ndx')
    cmd = 'echo -e "{index_Protein}|{index_LIG}\nq\n" | {gmx} make_ndx -n {ndxfile} -o {out_ndxfile}'.format(index_Protein=str(index_Protein), index_LIG=str(index_LIG),
                    gmx=gmx_bin, ndxfile=ndxfile, out_ndxfile=out_ndxfile)
    run_cmd(cmd)

    # if successful, the above command should make a new atomgroup called Protein_LIG,
    # ... for which we need to find the index
    index_Protein_LIG = get_ndx_group_index(out_ndxfile, 'Protein_LIG')
    print('index_Protein_LIG', index_Protein_LIG)

    ### Next let's build a Protein_LIG.gro from the setup npt.gro
    grofile = os.path.join(setup_rundir, 'npt.gro')
    out_grofile = os.path.join(setup_rundir, 'Protein_LIG.gro')
    cmd = 'echo -e "{index_Protein_LIG}\n" | {gmx} editconf -f {grofile} -n {out_ndxfile} -o {out_grofile}'.format(index_Protein_LIG=str(index_Protein_LIG),
                                 gmx=gmx_bin, out_ndxfile=out_ndxfile, grofile=grofile, out_grofile=out_grofile)
    run_cmd(cmd)

    ### Next, we create a topology file that only contains Protein and LIG
    """
    [ molecules ]
; Compound       #mols
LIG                  1
system1              1
HOH              16888
NA                  33
CL                  31
    """
    # find the topol.top topology file
    topfile = os.path.join(setup_rundir, 'topol.top')
    if not os.path.exists(topfile):
        print("Can't find topfile", topfile)
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), topfile)

    fin = open(topfile, 'r')
    toplines = fin.readlines()
    fin.close()

    while (toplines[-1].count('CL') + toplines[-1].count('NA') + toplines[-1].count('HOH') ) > 0:
        toplines.pop()

    out_topfile = topfile = os.path.join(setup_rundir, 'Protein_LIG.top')
    fout = open(topfile, 'w')
    fout.writelines(toplines)
    fout.close()
    print('Wrote:', out_topfile)


    #### Then, we write a quick-and-dirty *.mdp file, just so we can grompp to a tpr...
    mdp_txt = """;
; VARIOUS PREPROCESSING OPTIONS
; Preprocessor information: use cpp syntax.
; e.g.: -I/home/joe/doe -I/home/mary/roe
include                  = -I../top
; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
define                   = 
; RUN CONTROL PARAMETERS
integrator               = steep
emtol                    = 100.0
emstep                   = 0.03 ; was 0.05 prev to 3/20/2020
nsteps                   = 1000 ; 100000

comm-mode                = Linear
nstcomm                  = 1

; Output control
nstlog                   = 500		; every 1 ps
nstcalcenergy            = 1
nstenergy                = 50000        ; save edr every 100 ps
nstxout-compressed       = 50000	; save xtc coordinates every 100 ps
nstxout		 	 = 500000	; save coordinates every 1 ns
nstvout			 = 500000	; save velocities every 1 ns
compressed-x-precision	 = 1000
; 
; This selects the subset of atoms for the .xtc file. You can
; select multiple groups. By default all atoms will be written.
compressed-x-grps        = System

; Selection of energy groups
energygrps               = System

; Neighborsearching and short-range nonbonded interactions
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.0

; Electrostatics
cutoff-scheme            = verlet
coulombtype              = PME
coulomb-modifier         = none
rcoulomb                 = 1.0

; van der Waals
vdw-type                 = Cut-off
;;; vdw-modifier             = Potential-switch
rvdw-switch              = 0.9
rvdw                     = 1.0

; Apply long range dispersion corrections for Energy and Pressure 
DispCorr                 = EnerPres

; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
;fourier-nx               = 48
;fourier-ny               = 48
;fourier-nz               = 48
; EWALD/PME/PPPM parameters = 
pme_order                = 4
ewald_rtol               = 1e-05
ewald_geometry           = 3d
epsilon_surface          = 0

; Temperature coupling
tcoupl                   = v-rescale
nsttcouple               = 1
tc_grps                  = System
tau_t                    = 0.5
ref_t                    = 298.15
; Pressure coupling is on for NPT
pcoupl                   = no

; velocity generation
gen_vel                  = yes
gen-temp                 = 298.15
gen-seed                 = 5420 ; need to randomize the seed each time.

; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 4   ;12
lincs-iter               = 1   ;2
    """

    mdpfile = os.path.join(setup_rundir, 'Protein_LIG.mdp')
    fout = open(mdpfile, 'w')
    fout.write(mdp_txt)
    fout.close()


    # Grompp this topfile to get a *.tpr -- needed for Chimera; PDB transformations etc!
    out_tprfile = os.path.join(setup_rundir, 'Protein_LIG.tpr')
    cmd = '{gmx} grompp -f {mdpfile} -c {out_grofile} -p {out_topfile} -o {out_tprfile}'.format(
                    gmx=gmx_bin, mdpfile=mdpfile, out_grofile=out_grofile, out_topfile=out_topfile, out_tprfile=out_tprfile)
    run_cmd(cmd)
    print('Wrote:', out_tprfile)


  

def get_ndx_group_index(ndxfile, index_strings):
    """Gets the atom group index for the given strings.
    
    INPUT
    ndxfile        - index file *.ndx filename
    index_strings  - either a single string, or a *list* of strings
                     This string MUST be in the form 'Example' 
                     where the field in the *.ndx file is '[ Example ]'
    
    RETURNS
    results        - if a single string supplied, the index (int)
                   - elif a list of strings supples, a list of ints

    """

    assert ((type(index_strings) == str) | (type(index_strings) == list))

    # get all the lines with index headers
    fin = open(ndxfile, 'r')
    index_lines = [ line.strip() for line in fin.readlines() if line.count('[') > 0 ]
    fin.close()

    print('index_strings', index_strings)
    print('index_lines', index_lines)

    if type(index_strings) == str:
        return index_lines.index( '[ '+index_strings+' ]' )
    if type(index_strings) == list:
        return [ index_lines.index( '[ '+s+' ]' ) for s in index_strings]

    return None

"""
gmx trjcat -o all.xtc -f results?/traj_comp.xtc  results??/traj_comp.xtc -cat
gmx trjconv -s ../../../PROJ14366-RUN762-setup/Protein_LIG.tpr -o all_cluster.xtc -f all.xtc -pbc cluster

gmx trjcat -o all.xtc -f results?/traj_comp.xtc  results??/traj_comp.xtc -cat
gmx trjconv -s ../../../PROJ14366-RUN762-setup/Protein_LIG.tpr -o all_cluster.xtc -f all.xtc -pbc cluster
"""


if __name__ == "__main__":
    #import doctest
    #doctest.testmod()

    build_Protein_LIG_tpr('/home/server/server2/projects/p14399/RUN0')
