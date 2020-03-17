import os, sys, glob
import numpy as np


usage = """
    Usage:   make_fep_ready.py [input rundir] [output rundir] ligonly

    use "ligonly" keyword if you are readying ligand-only simulatons

NOTES
    * This should be used with an installation of GROMACS 5.0.4, as on the FAH servers!!!!
    * It assumes you have a structure `conf.gro` that represents the whole system,
      a topology `topol.top`

EXAMPLE
    $ python make_fep_ready.py ../100_ligands/RUN1 ./RUN1
or 
    $ python make_fep_ready.py ../100_ligands_noreceptor/RUN1 ./RUN1 ligonly
"""

if len(sys.argv) < 3:
    print(usage)
    sys.exit(1)

# parse the input arguments
in_rundir = sys.argv[1]
out_rundir = sys.argv[2]
if not os.path.exists(out_rundir):
    os.mkdir(out_rundir)
ligand_only = False
if len(sys.argv) > 3:
    if sys.argv[3].count('only') > 0:
        ligand_only = True


############# convert the grofile ###############3
# Find the conf.gro file
in_grofile = os.path.join(in_rundir, 'conf.gro')
if not os.path.exists(in_grofile):
    print("Can't find grofile", in_grofile, '! Exiting.')
else:
    print('Found:', in_grofile)

out_grofile = os.path.join(out_rundir, 'minimize_me.gro')

# read in the lines
fin = open(in_grofile, 'r')
gro_contents = fin.read()
fin.close()

print('Writing to', out_grofile, '...')
fout = open(out_grofile, 'w')
gro_contents = gro_contents.replace('UNL', 'LIG')
fout.write( gro_contents )
fout.close()
print('...Done.')

############### convert to topfile ###############
in_topfile = os.path.join(in_rundir, 'topol.top')
if not os.path.exists(in_topfile):
    print("Can't find topfile", in_topfile, '! Exiting.')
else:
    print('Found:', in_topfile)

out_topfile = os.path.join(out_rundir, 'topol.top')

# read in the lines
fin = open(in_topfile, 'r')
top_contents = fin.read()
fin.close()

print('Writing to', out_topfile, '...')
fout = open(out_topfile, 'w')
top_contents = top_contents.replace('UNL', 'LIG')
fout.write( top_contents )
fout.close()
print('...Done.')


############## make an index file ###########

# Make the "default" index file that gmx makes from a grofile
ndxfile = os.path.join(out_rundir, 'index.ndx')

## if it already exists, back it up
if os.path.exists(ndxfile):
    print('Index file exists! Backing up... ', end='')
    bkup_cmd = 'mv {ndxfile} {ndxfile}.bkup'.format(ndxfile=ndxfile) 
    print(bkup_cmd) 
    os.system(bkup_cmd)

GMX_BIN = os.environ['GMXBIN']
cmd = 'echo "q\\n" | {GMX_BIN}/gmx make_ndx -f {out_grofile} -o {ndxfile}'.format(GMX_BIN=GMX_BIN, out_grofile=out_grofile, ndxfile=ndxfile)
os.system(cmd)

if not ligand_only:

    ## Finally, we have to select two atoms for the restraint 

    ### Pick '  165MET     CA' as the atom on the protein
    gro_lines = [ line for line in gro_contents.split('\n') if line.count('  165MET     CA') > 0 ]
    print('gro_lines', gro_lines)
    if len(gro_lines)  == 0:
        print("Can't find 165MET     CA in the grofile -- what gives????  Exiting.")
        sys.exit(1)
    protein_atomnum = int( (gro_lines[0].strip()).split()[2] )
    protein_xyz = np.array( [float(s) for s in (gro_lines[0].strip()).split()[3:]] )   # units nm
    print('protein_atomnum', protein_atomnum)
    print('protein_xyz', protein_xyz)
    
    ### Pick the closest carbon in the LIG to the protein_atom
    lig_lines = [ line for line in gro_contents.split('\n') if ( (line.count('LIG') > 0) and (line.count(' C')>0) ) ]
    print('lig_lines', lig_lines)
    if len(lig_lines)  == 0:
        print("Can't find LIG and C in grofile -- what gives????  Exiting.")
        sys.exit(1)
    ligand_positions = []
    for i in range(len(lig_lines)):
        ligand_positions.append( [float(s) for s in (lig_lines[i].strip()).split()[3:]] )
    ligand_positions = np.array(ligand_positions)
    print('ligand_positions', ligand_positions)
    ligand_vectors = ligand_positions - np.tile(protein_xyz, (ligand_positions.shape[0],1))  # N, 3
    print('ligand_vectors', ligand_vectors)
    ligand_distances = ligand_vectors[:,0]**2 + ligand_vectors[:,1]**2 +ligand_vectors[:,2]**2  
    print('ligand_distances', ligand_distances)
    i_closest = np.argmin(ligand_distances)
    print('i_closest', i_closest)
    print('lig_lines[i_closest]', lig_lines[i_closest])
    ligand_atomnum = int( (lig_lines[i_closest].strip()).split()[2] )
    print('ligand_atomnum', ligand_atomnum)
   
    fout = open(ndxfile, 'a')  # append these new lines to the index file...
    index_txt = """[ a1-Protein ]
      {protein_atomnum} 
    [ a2-Ligand ]
      {ligand_atomnum}
    [ Restraint-Distance ]
      {protein_atomnum} {ligand_atomnum}
    """.format(protein_atomnum=protein_atomnum, ligand_atomnum=ligand_atomnum)
    fout.write(index_txt)
    fout.close()


### do a very quick minimization to make sure we're all set to run

min_mdpfile = os.path.join(out_rundir, 'quickmin.mdp')
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
emstep                   = 0.05
nsteps                   = 100000 ; 100000

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
compressed-x-grps        = Protein LIG

; Selection of energy groups
energygrps               = Protein non-Protein
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

; pulling parameters
pull                     = umbrella
pull_ngroups             = 2
pull_ncoords             = 1
pull_group1_name         = a1-Protein
pull_group2_name         = a2-Ligand
pull-geometry            = direction-periodic
pull_coord1_groups       = 1 2
pull-dim                 = Y Y Y
pull_coord1_rate         = 0.00
pull_coord1_k            = 200.0
pull-start               = yes   ;    no
pull-coord1-init         = 0.375
pull-nstxout             = 500   ; 1 ps
pull-nstfout             = 500   ; 1 ps



"""
fout = open(min_mdpfile, 'w')
fout.write(mdpfile_text)
fout.close()

# run the minimization
os.system( 'gmx grompp -c {out_grofile} -f {out_rundir}/quickmin.mdp -p {out_topfile} -n {ndxfile} -o {out_rundir}/min.tpr'.format(out_grofile=out_grofile, out_rundir=out_rundir, out_topfile=out_topfile, ndxfile=ndxfile) )
os.system( 'gmx mdrun -v -s {out_rundir}/min.tpr -c {out_rundir}/npt.gro'.format(out_rundir=out_rundir) )

# print some commands to test the script
testing_cmds = "### To test this project, try: ###\n"

if ligand_only:
    testing_cmds += "python ../scripts/create_ee_mdp.py {out_rundir}/npt.gro {out_topfile} {ndxfile} {out_rundir}/prod.mdp ligonly\n".format(out_topfile=out_topfile, ndxfile=ndxfile, out_rundir=out_rundir)

    testing_cmds += """mkdir test
gmx grompp -c {out_rundir}/npt.gro -f {out_rundir}/prod.mdp -p {out_topfile} -n {ndxfile} -o test/testme.tpr -po mdout.mdp -maxwarn 1
cd test
gmx mdrun -nt 1 -v -s testme.tpr""".format(out_topfile=out_topfile, ndxfile=ndxfile, out_rundir=out_rundir)

else:
    testing_cmds += "python ../scripts/create_ee_mdp.py {out_rundir}/npt.gro {out_topfile} {ndxfile} {out_rundir}/prod.mdp\n".format(out_topfile=out_topfile, ndxfile=ndxfile, out_rundir=out_rundir)

    testing_cmds += """mkdir test
gmx grompp -c {out_rundir}/npt.gro -f {out_rundir}/prod.mdp -p {out_topfile} -n {ndxfile} -o test/testme.tpr -po mdout.mdp -maxwarn 1
cd test
gmx mdrun -v -s testme.tpr""".format(out_topfile=out_topfile, ndxfile=ndxfile, out_rundir=out_rundir) 

print(testing_cmds)