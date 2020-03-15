import os, sys, glob
import numpy as np


usage = """
    Usage:   make_fep_ready.py [input rundir] [output rundir]

NOTES
    * This should be used with an installation of GROMACS 5.0.4, as on the FAH servers!!!!
    * It assumes you have a structure `conf.gro` that represents the whole system,
      a topology `topol.top`

EXAMPLE
    $ python make_fep_ready.py ../100_ligands/RUN1 ./RUN1

"""

if len(sys.argv) < 3:
    print(usage)
    sys.exit(1)

# parse the input arguments
in_rundir = sys.argv[1]
out_rundir = sys.argv[2]
if not os.path.exists(out_rundir):
    os.mkdir(out_rundir)


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
emstep                   = 0.01
nsteps                   = 500 ; 100000
cutoff-scheme            = Verlet
energygrps               = System
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 0.9
coulombtype              = PME
rcoulomb-switch          = 0
rcoulomb                 = 0.9
epsilon_r                = 1
epsilon_rf               = 1
vdw-type                 = Cut-off
rvdw-switch              = 0
rvdw                     = 0.9
fourierspacing           = 0.12
pme_order                = 4
ewald_rtol               = 1e-5
ewald_geometry           = 3d
epsilon_surface          = 0
optimize_fft             = no
;tcoupl                   = V-rescale
compressibility          = 4.5e-5
ref_p                    = 1.0
constraints              = hbonds
constraint-algorithm     = Lincs
continuation             = no
Shake-SOR                = no
shake-tol                = 0.0001
lincs-order              = 4
lincs-iter               = 1
lincs-warnangle          = 30
morse                    = no
nwall                    = 0
wall_type                = 9-3
wall_r_linpot            = -1
wall_atomtype            = 
wall_density             = 
wall_ewald_zfac          = 3
"""
fout = open(min_mdpfile, 'w')
fout.write(mdpfile_text)
fout.close()

# run the minimization
os.system( 'gmx grompp -c {out_grofile} -f {out_rundir}/quickmin.mdp -p {out_topfile} -n {ndxfile} -o {out_rundir}/min.tpr'.format(out_grofile=out_grofile, out_rundir=out_rundir, out_topfile=out_topfile, ndxfile=ndxfile) )
os.system( 'gmx mdrun -v -s {out_rundir}/min.tpr -c {out_rundir}/npt.gro'.format(out_rundir=out_rundir) )

# print some commands to test the script
print()
print('### To test this project, try: ###')
print("""
python ../scripts/create_ee_mdp.py {out_rundir}/npt.gro {out_topfile} {ndxfile} {out_rundir}/prod.mdp
mkdir test
gmx grompp -c {out_rundir}/npt.gro -f {out_rundir}/prod.mdp -p {out_topfile} -n {ndxfile} -o test/testme.tpr -po mdout.mdp -maxwarn 1
cd test
gmx mdrun -v -s testme.tpr
""".format(out_grofile=out_grofile, out_topfile=out_topfile, ndxfile=ndxfile, out_rundir=out_rundir) )





