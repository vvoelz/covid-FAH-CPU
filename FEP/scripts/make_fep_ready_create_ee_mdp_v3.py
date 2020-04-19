import os, sys, glob, stat
import subprocess
import errno

import numpy as np

from expanded_v2 import *

###############################################################
### First thing's first: do we have a GROMACS install or what?
global GMX_BIN
try:
    GMX_BIN = os.environ['GMXBIN']
except:
    GMX_BIN = '/usr/local/gromacs/bin'
if not os.path.exists(GMX_BIN):
    print('Cannot find the GROMACS installation!  Please set environment variable GMX_BIN')
    sys.exit(1)
###############################################################


# Subroutines and Functions

def rename_grofile_residues(in_grofile, out_grofile, old_res_names=['UNK', 'UNL'])
    """This function will read in a grofile and change 'UNK', 'UNL', etc. residue names to 'LIG'. """

    ############# convert the grofile ###############3
    # Find the conf.gro file
    in_grofile = os.path.join(in_rundir, 'conf.gro')
    if not os.path.exists(in_grofile):
        print("Can't find grofile", in_grofile)
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), ndxfile)
    else:
        print('Found:', in_grofile)

    out_grofile = os.path.join(out_rundir, 'minimize_me.gro')

    # read in the lines
    fin = open(in_grofile, 'r')
    gro_contents = fin.read()
    fin.close()

    print('Writing to', out_grofile, '...')
    fout = open(out_grofile, 'w')
    for resname in old_res_names:
        gro_contents = gro_contents.replace(resname, 'LIG')
    fout.write( gro_contents )
    fout.close()
    print('...Done.')

    return

def modify_topfile(in_topfile, out_topfile):
    """modify the topfile"""

    if not os.path.exists(in_topfile):
        print("Can't find topfile", in_topfile, '! Exiting.')
    else:
        print('Found:', in_topfile)

    # read in the lines
    fin = open(in_topfile, 'r')
    top_lines = fin.readlines()
    fin.close()

    for i in range(len(top_lines)):

    # Change names from UNL to LIG
    top_lines[i] = top_lines[i].replace('UNL', 'LIG').replace('UNK', 'LIG')

    # Convert any hydrogen mass to 4.000 amu
    ### NOTE there is no error-checking here -- to do: make this more robust
    if top_lines[i].count('H') > 0:   # We ASSUME that atom names have an 'H'
        top_lines[i] = top_lines[i].replace('1.007947', '4.000000')

    print('Writing to', out_topfile, '...')
    fout = open(out_topfile, 'w')
    fout.writelines( top_lines )
    fout.close()
    print('...Done.')


def active_site_restraint_info(grofile, residues=['ALA', 'VAL'], cutoff=1.0, verbose=True):
    """Get information about protein and ligand atom groups to restrain.
    
    INPUT
    grofile        - an input grofile with the LIG residue correctly named
    
    PARAMETERS
    residues       - list of three-letter amino acid residue codes to search. (Default: ['ALA', 'VAL'])
    cutoff         - the cutoff distance (in nm) threshold to include residue alphha-carbons in the restraint list.
                     The distance is computed as d(x2 - x1) where:
                     
                         x1 is the center of mass (COM) of all carbons in the ligand,
                         x2 is each alpha-carbon in the specified residues.

                     (Default: 1.0 nm)
    RETURNS
    protein_indices  - a list of protein CA indices (in the gmx convention, starting at 1,  gmx) for atom group 1
    ligand_indices   - a list of ligand C indices ((in the gmx convention, starting at 1,  gmx) for atom group 2
    com_distance   - the distance between the center of masses of the protein_atoms and ligand_atoms
    """

    # read in the grofilelines
    fin = open(grofile, 'r')
    gro_contents = fin.read()
    fin.close()

    ##########################################################
    # First, let's find all the carbon atoms in the ligand

    ### get all the LIG lines with carbon
    lig_lines = [ line for line in gro_contents.split('\n') if ( (line.count('LIG') > 0) and (line.count(' C')>0) ) ]
    print('lig_lines', lig_lines)
    if len(lig_lines)  == 0:
        print("Can't find LIG and C in grofile -- what gives????  Exiting.")
        sys.exit(1)

    ### compile the ligand C indices 
    ligand_indices = []
    for i in range(len(lig_lines)):
        ligand_indices.append( [int(s) for s in (lig_lines[i].strip()).split()[2]] )
    print('ligand_indices', ligand_indices)

    ### compile the ligand C positions
    ligand_positions = []
    for i in range(len(lig_lines)):
        ligand_positions.append( [float(s) for s in (lig_lines[i].strip()).split()[3:6]] )
    ligand_positions = np.array(ligand_positions)
    print('ligand_positions', ligand_positions)

    ### compute the COM of the ligand C atoms
    com_ligand = np.mean(ligand_positions, axis=0)
    print('com_ligand', com_ligand)

    ##########################################################
    # Next, let's find all of the residue alpha-carbon atom indices that are within the cutoff

    protein_lines = []
    for res in residues:
        protein_lines += [ line for line in gro_contents.split('\n') if line.count('{res}     CA'.format(res=res)) > 0 ]
    if len(protein_lines) == 0:
        print("Can't find any ALA of VAL residues????   - what is this crazy protein strutcure?  formatting problems?  Exiting.")
        sys.exit(1)

    protein_indices = []
    protein_positions = []
    for protein_line in protein_lines:
        protein_index = int( (protein_line.strip()).split()[2] )
        protein_xyz = np.array( [float(s) for s in (protein_line.strip()).split()[3:6]] )   # units nm

        # is the protein_xyz closer than the cutoff to the com_ligand?
        vec = (protein_xyz - com_ligand)
        distance = ((vec**2).sum())**0.5

        if distance <= cutoff:
            protein_indices.append(protein_index)
            protein_positions.append(protein_xyz)

    protein_positions = np.array(protein_positions)   # convert to an array!

    print('protein_indices', protein_indices)
    print('protein_positions', protein_positions)

    # Find the center of mass of the selected protein_positions
    com_protein = protein_positions.mean(axis=0)

    # compute the com_distance between the two groups
    vec = (com_protein - com_ligand)
    com_distance =  ((vec**2).sum())**0.5

    return protein_indices, ligand_indices, com_distance



def quick_minimize(input_grofile, topfile, ndxfile, output_grofile, workdir, babysteps=False):
    """Perform a quick minimization.

    A 'quickmin.mdp' and 'quickmin.tpr' will be written to the workdir, and
    the final minimized structure will be written to out_grofile.

    OPTIONS
    babysteps    - if True, use a smaller emstep in the minimization (Default: False)."""

    min_mdpfile = os.path.join(workdir, 'quickmin.mdp')
    min_tprfile = os.path.join(workdir, 'quickmin.tpr')

    if babysteps:
        emstep = 0.003
    else:
        emstep = 0.03

    mdpfile_text = """;       generated from expanded_ensemble_mdpfile() on Sat Mar 14 21:04:11 2020
;
;
; VARIOUS PREPROCESSING OPTIONS
; Preprocessor information: use cpp syntax.
; e.g.: -I/home/joe/doe -I/home/mary/roe
include                  = -I../top
; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
define                   =
; RUN CONTROL PARAMETERS
integrator               = steep
emtol                    = 100.0
emstep                   = {emstep} ; was 0.05 prev to 3/20/2020
nsteps                   = 100 ; Trying this Apr 3, 2020 VAV  100000

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
compressed-x-grps        = LIG

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

""".format(emstep=emstep)

    fout = open(min_mdpfile, 'w')
    fout.write(mdpfile_text)
    fout.close()

    # run the minimization
    os.system( '{GMX_BIN}/gmx grompp -c {input_grofile} -f {min_mdpfile} -p {topfile} -n {ndxfile} -o {min_tprfile}'.format(GMX_BIN=GMX_BIN,
                input_grofile=input_grofile, min_mdpfile=min_mdpfile, topfile=topfile, ndxfile=ndxfile, min_tprfile=min_tprfile) )
    os.system( '{GMX_BIN}/gmx mdrun -v -s {min_tprfile} -c {output_grofile}'.format(GMX_BIN=GMX_BIN,
                min_tprfile=min_tprfile, output_grofile=output_grofile) )

    # print some commands to test the script
    testing_cmds = "### To test this project, try: ###\n"

    if ligand_only:
        testing_cmds += "python ../scripts/create_ee_mdp.py {out_rundir}/npt.gro {out_topfile} {ndxfile} {out_rundir}/prod.mdp ligonly\n".format(out_topfile=out_topfile, ndxfile=ndxfile, out_rundir=out_rundir)

        testing_cmds += """mkdir test
{GMX_BIN}/gmx grompp -c {out_rundir}/npt.gro -f {out_rundir}/prod.mdp -p {out_topfile} -n {ndxfile} -o test/testme.tpr -po mdout.mdp -maxwarn 1
cd test
{GMX_BIN}/gmx mdrun -nt 1 -v -s testme.tpr""".format(GMX_BIN=GMX_BIN, out_topfile=out_topfile, ndxfile=ndxfile, out_rundir=out_rundir)

    else:
        testing_cmds += "python ../scripts/create_ee_mdp.py {out_rundir}/npt.gro {out_topfile} {ndxfile} {out_rundir}/prod.mdp\n".format(out_topfile=out_topfile, ndxfile=ndxfile, out_rundir=out_rundir)

        testing_cmds += """mkdir test
{GMX_BIN}/gmx grompp -c {out_rundir}/npt.gro -f {out_rundir}/prod.mdp -p {out_topfile} -n {ndxfile} -o test/testme.tpr -po mdout.mdp -maxwarn 1
cd test
{GMX_BIN}/gmx mdrun -v -s testme.tpr""".format(GMX_BIN=GMX_BIN, out_topfile=out_topfile, ndxfile=ndxfile, out_rundir=out_rundir)

    print(testing_cmds)



def create_ndxfile(grofile, ndxfile, verbose=True):
    """Create an ndxfile from the given grofile.  If it already exists, make a backup."""

    ## if it already exists, back it up
    if os.path.exists(ndxfile):
        print('Index file exists! Backing up... ', end='')
        bkup_cmd = 'mv {ndxfile} {ndxfile}.bkup'.format(ndxfile=ndxfile)
        print(bkup_cmd)
        os.system(bkup_cmd)

    cmd = 'echo "q\\n" | {GMX_BIN}/gmx make_ndx -f {grofile} -o {ndxfile}'.format(GMX_BIN=GMX_BIN, grofile=grofile, ndxfile=ndxfile)
    os.system(cmd)


#############################################

# Main

if __name__ == '__main__':

    print('Hello world!')
    sys.exit(1)

    import argparse, textwrap

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\

    Process COVID-19 protease drug screening setup directories to make FEP-ready

    THIS is version v3 -- 
    The protocol is to find all ALA and VAL alpha-carbons  within 10.0 Angstroms from the ligand and use those groups

    NOTES
    * This should be used with an installation of GROMACS 5.0.4, as on the FAH servers!!!!
    * It assumes you have a structure `conf.gro` that represents the whole system,
      and a topology `topol.top`

      use "--ligonly" ONLY if you are readying ligand-only simulatons

    EXAMPLE
    $ python make_fep_ready_create_ee_mdp_v3.py ../100_ligands/RUN1 ./RUN1
or
    $ python make_fep_ready_create_ee_mdp_v3.py ../100_ligands_noreceptor/RUN1 ./RUN1 --ligonly

      ''' ))

    parser.add_argument('in_rundir', type=str, help='the input rundir')
    parser.add_argument('out_rundir', type=str, help='the output rundir')
    parser.add_argument('--ligonly', dest='ligand_only', action='store_true',
                    help='Specify this is a ligand-only simulation')
    parser.add_argument('--babysteps', dest='babysteps', action='store_true',
                    help='Use "baby steps" in the energy minimization to avoid errors')
    parser.add_argument('--moonshot', dest='moonshot', action='store_true',
                    help='Specify the structure comes from a moonshot simulation.  The resnum in conf.gro is off by +2!')
    parser.add_argument('--agg', dest='agg', action='store_true',
                    help='Specify the structure comes from an AGG simulation.  The resnum in conf.gro is off by +1!')

    args = parser.parse_args()
    print('args.in_rundir', args.in_rundir)
    print('args.out_rundir', args.out_rundir)
    print('args.ligand_only', args.ligand_only)
    print('args.babysteps', args.babysteps)


    # parse the input arguments
    in_rundir = args.in_rundir
    out_rundir = args.out_rundir
    if not os.path.exists(out_rundir):
        os.mkdir(out_rundir)
    ligand_only = args.ligand_only
    babysteps = args.babysteps
    moonshot = args.moonshot
    agg = args.agg


    ### convert the grofile to FEP-ready by renaming the ligand residues to LIG ###
    in_grofile = os.path.join(in_rundir, 'conf.gro')
    out_grofile = os.path.join(out_rundir, 'minimize_me.gro')
    rename_grofile_residues(in_grofile, out_grofile)


    ### convert the topfile to FEP-ready by changing resnames and hydrogen masses ###
    in_topfile = os.path.join(in_rundir, 'topol.top')
    out_topfile = os.path.join(out_rundir, 'topol.top')
    modify_topfile(in_topfile, out_topfile)


    ############## make a default index file ###########
    ndxfile = os.path.join(out_rundir, 'index.ndx')
    create_ndxfile(out_grofile, ndxfile)

    ###  NEW in v3 -- let's find a set of protein and ligand indices whose COM distance we will restrain
    if not ligand_only:
        protein_indices, ligand_indices, com_distance = active_site_restraint_info(out_grofile, residues=['ALA', 'VAL'], cutoff=0.10, ligand_atoms='carbon')

    ### add these pull groups to the index file
    pull_group_1 = ' '.join( [str(i) for i in protein_indices ])
    pull_group_2 = ' '.join( [str(i) for i in ligand_indices ])
    fout = open(ndxfile, 'a')  # append these new lines to the index file...
    index_txt = """[ a1-Protein ]
  {pull_group_1} 
[ a2-Ligand ]
  {pull_group_2}
    """.format(pull_group_1=pull_group_1, pull_group_2=pull_group_2)
    fout.write(index_txt)
    fout.close()

    ### do a very quick minimization to make sure we're all set to run
    minimized_grofile = os.path.join(out_rundir, 'npt.gro')
    quick_minimize(out_grofile, out_topfile, ndxfile, minimized_grofile, out_rundir)

    ### Create an expanded-ensemble mdpfile the correct pull groups and equilibrium distance 
    if not ligand_only:

        # To tether the ligand to the protein, we need to have prepared an index file with
        # the following atom groups, for example:
        """
        [ a1-Protein ]
        678 690 .....
        [ a2-Ligand ]
        1564 1565 1566 ...
        """
        # NOTE that for our tool chain to work smoothly, we should always use these DEFAULT
        # atom group names

        # use the expanded_ensemble_mdpfile() class to create an initial mdpfile
        e = expanded_ensemble_mdpfile(ligand_only = ligand_only,
                                      pull_group1_name = 'a1-Protein',
                                      pull_group2_name = 'a2-Ligand',
                                      pull_coord1_init  = com_distance)
    else:
        e = expanded_ensemble_mdpfile(ligand_only = ligand_only)
    
    e.write_to_filename(mdpfile)



