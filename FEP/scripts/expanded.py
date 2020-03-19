import os, sys
import time
import numpy as np

# Classes and Functions for preparing and continuing expanded-ensemble simulations

class expanded_ensemble_mdpfile(object):
    """A class stucture for creating and writing standardized *.mdp files to do
    expanded-ensemble GROMACS 5.0.4 simualtions ."""

    def __init__(self,  ligand_only=False,
                        couple_moltype = 'LIG',
			pull_group1_name = 'a1-Protein', 
                        pull_group2_name = 'a2-Ligand', 
                        pull_coord1_k    = 50.0, 
                        pull_coord1_init  = 0.4,   # an initial guess; gets measured in make_fep_ready.py
                        fep_lambdas      = np.array( 40*[0.0] ),
                        coul_lambdas     = np.array( [ 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                                                       0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
                                                       1.0, 1.0,  1.0, 1.0,  1.0, 1.0,  1.0, 1.0,  1.0, 1.0, 
                                                       1.0, 1.0,  1.0, 1.0,  1.0, 1.0,  1.0, 1.0,  1.0, 1.0 ] ), 
                        vdw_lambdas      = np.array( [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                                       0.00, 0.10, 0.20, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.63,
                                                       0.66, 0.69, 0.72, 0.75, 0.78, 0.81, 0.84, 0.88, 0.92, 1.00 ] ), 
                        init_lambda_weights = np.zeros(40),
                        init_lambda_state   = 0,
                        wl_increment_in_kT  = 10.0 ):



        """Initialize the class with default options.

        OPTIONS

        ligand_only -   If True, then the mdp will be written WITHOUT the pull code and without
                        reference to any protein.

        [ to do -- add other options here ] 

        NOTES

        The *.mdp files that this class reads and write MUST have standard naming conventions

        For the ligand,
            * the Gromacs topology `*.top` MUST have a moleculetype `LIG`
            * the index file `*.ndx` MUST have an index group `LIG` 

        ### Protein-ligand simulations ###

        For a protein-ligand simulation, there is a tether between the ligand to the protein,
        and we need to have prepared an index file with the atom groups `a1-Protein`, `a2-Ligand` and
        `Restraint-Distance`, define as in the following example:
            ....
            [ a1-Protein ]
            678
            [ a2-Ligand ]
            1564
            [ Restraint-Distance ]
            678 1564

        Moreover, there needs to be atom groups `Protein` and `non-Protein`, which the thermostat will control separately.

        ### Ligand-only simulations ###

        There needs to be atom groups `Water` and `non-Water`, which the thermostat will control separately.


        OPTIONS

        ligand_only -   If True, then the mdp will be written WITHOUT the pull code and without
                        reference to any protein."""

        self.ligand_only      = ligand_only 
        self.couple_moltype   = couple_moltype 
        self.pull_group1_name = pull_group1_name
        self.pull_group2_name = pull_group2_name
        self.pull_coord1_k    = pull_coord1_k     # the spring constant in kJ/nm^2
        self.pull_coord1_init = pull_coord1_init  # the equilibrium distance between the pull group 1 and 2
        self.fep_lambdas      = fep_lambdas
        self.nlambdas         = len(self.fep_lambdas)
        self.fep_lambdas_string   = ' '.join(['%2.3f'%fep_lambdas[i] for i in range(self.nlambdas)])

        self.init_lambda_weights  = init_lambda_weights
        self.init_lambda_weights_string =  ' '.join(['%2.3f'%init_lambda_weights[i] for i in range(self.nlambdas)])

        self.coul_lambdas     = coul_lambdas
        self.coul_lambdas_string   = ' '.join(['%2.3f'%coul_lambdas[i] for i in range(self.nlambdas)])

        self.vdw_lambdas     = vdw_lambdas
        self.vdw_lambdas_string   = ' '.join(['%2.3f'%vdw_lambdas[i] for i in range(self.nlambdas)])

        self.init_lambda_state    = init_lambda_state
        self.wl_increment_in_kT   = wl_increment_in_kT

        # Here's a random number seed, in case we need it:
        self.randseed         = np.random.randint(1e5)


    def report(self):
        """Print a report of the settings."""

        print('self.ligand_only', self.ligand_only)
        print('self.couple_moltype', self.couple_moltype)
        print('self.pull_group1_name', self.pull_group1_name)
        print('self.pull_group2_name', self.pull_group2_name)
        print('self.pull_coord1_k', self.pull_coord1_k)
        print('self.pull_coord1_init', self.pull_coord1_init)
        print('self.fep_lambdas', self.fep_lambdas)
        print('self.nlambdas', self.nlambdas)
        print('self.fep_lambdas_string', self.fep_lambdas_string)
        print('self.init_lambda_weights', self.init_lambda_weights)
        print('self.init_lambda_weights_string', self.init_lambda_weights_string)
        print('self.init_lambda_state', self.init_lambda_state)
        print('self.wl_increment_in_kT', self.wl_increment_in_kT)



    def read_parms_from_mdpfile(self, input_mdpfile, VERBOSE=True):
        """Read in the key parameters from an existing mdpfile.

        WARNING -- there's no error checking here.   Something for the future..."""

        # Read in all lines of the mdpfile
        fin = open(input_mdpfile, 'r')
        lines = fin.readlines() 
        fin.close()

        # Search for keywords
        print('---- Searching for keywords ----')
        for line in lines:

            # REMOVE COMMENTS -- anything before ';'
            line = line.split(';')[0]

            # Use the temperature control groups to detect ligand_only
            if line.count('tc-grps') > 0:
                if line.count('Protein'):
                    self.ligand_only = False
                else:
                    self.ligand_only = True
                print('Setting ligand_only -->', self.ligand_only)

            # couple-moltype
            if line.count('couple-moltype') > 0:   # couple-moltype           = {couple_moltype}
                fields = line.split('=') 
                self.couple_moltype = fields[1].strip()
                print('Setting self.couple_moltype -->', self.couple_moltype)

            # pull_group1_name
            if line.count('pull_group1_name') > 0:   # pull_group1_name         = {pull_group1_name}
                fields = line.split('=')
                self.pull_group1_name = fields[1].strip()
                print('Setting self.pull_group1_name -->', self.pull_group1_name)

            # pull_group2_name
            if line.count('pull_group2_name') > 0:   # pull_group2_name         = {pull_group2_name}
                fields = line.split('=')
                self.pull_group2_name = fields[1].strip()
                print('Setting self.pull_group2_name -->', self.pull_group2_name)

            # pull_coord1_k
            if line.count('pull_coord1_k') > 0:   # pull_coord1_k            = {pull_coord1_k}
                fields = line.split('=')
                self.pull_coord1_k = fields[1].strip()
                print('Setting self.pull_coord1_k -->', self.pull_coord1_k)

            # fep_lambdas
            if line.count('fep-lambdas') > 0:   # fep-lambdas              = {fep_lambdas_string}
                fields = line.split('=')
                self.fep_lambdas_string = fields[1].strip()
                print('Setting self.fep_lambdas_string -->', self.fep_lambdas_string)
                self.fep_lambdas = np.array( [ float(s) for s in self.fep_lambdas_string.split() ] )
                print('Setting self.fep_lambdas --> ', self.fep_lambdas)
                self.nlambdas = self.fep_lambdas.shape[0]
                print('Setting self.nlambdas --> ', self.nlambdas)

            # init_lambda_state
            if line.count('init-lambda-state') > 0:  # init-lambda-state        = {init_lambda_state}
                fields = line.split('=')
                self.init_lambda_state = int(fields[1].strip())
                print('Setting self.init_lambda_state -->', self.init_lambda_state)

            # init_lambda_weights
            if line.count('init-lambda-weights') > 0:  # init-lambda-weights      = {init_lambda_weights}
                fields = line.split('=')
                self.init_lambda_weights_string = fields[1].strip()
                print('Setting self.init_lambda_weights_string -->', self.init_lambda_weights_string)
                self.init_lambda_weights = np.array( [ float(s) for s in self.init_lambda_weights_string.split() ] )
                print('Setting self.init_lambda_weights -->', self.init_lambda_weights)

            # wl_increment_in_kT
            if line.count('init-wl-delta') > 0:  # init-wl-delta            = {wl_increment_in_kT}  ; kT -- default is 1.0
                fields = line.split('=')
                self.wl_increment_in_kT = float(fields[1].strip())
                print('Setting self.wl_increment_in_kT -->', self.wl_increment_in_kT)


        # print a final report
        print('\nThe rest of the settings have their defaults retained:')
        self.report()

         

    def mdpfile_text(self):
        """Builds and returns a string corresponding to the mdpfile contents."""

        self.header_desc = ";       generated from expanded_ensemble_mdpfile() on %s \n;\n;\n"%time.asctime()

        self.mdp_text = ''

        self.mdp_text +=  """; Run control
integrator               = md-vv
tinit                    = 0
dt                       = 0.004
nsteps                   = 250000       ; 1 ns
comm-mode                = Linear
nstcomm                  = 1

; Output control
nstlog                   = 250		; every 1 ps
nstcalcenergy            = 1
nstenergy                = 25000        ; save edr every 100 ps
nstxout-compressed       = 25000	; save xtc coordinates every 100 ps
nstxout		 	 = 250000	; save coordinates every 1 ns
nstvout			 = 250000	; save velocities every 1 ns
compressed-x-precision	 = 1000
"""


        if self.ligand_only:
            self.mdp_text +=  """; 
; This selects the subset of atoms for the .xtc file. You can
; select multiple groups. By default all atoms will be written.
compressed-x-grps        = LIG

; Selection of energy groups
energygrps               = Water non-Water
"""
        else:
            self.mdp_text +=  """; 
; This selects the subset of atoms for the .xtc file. You can
; select multiple groups. By default all atoms will be written.
compressed-x-grps        = Protein LIG

; Selection of energy groups
energygrps               = Protein non-Protein
"""


        self.mdp_text +=  """; Neighborsearching and short-range nonbonded interactions
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 0.9

; Electrostatics
cutoff-scheme            = verlet
coulombtype              = PME
rcoulomb                 = 0.9

; van der Waals
vdw-type                 = Cut-off
vdw-modifier             = Potential-switch
rvdw-switch              = 0.89      ;    0.9
rvdw                     = 0.9

; Apply long range dispersion corrections for Energy and Pressure 
; NO -- we're doing NVT
; DispCorr                 = EnerPres

fourierspacing           = 0.10
pme_order                = 4
ewald_rtol               = 1e-6
ewald_geometry           = 3d
epsilon_surface          = 0


; Temperature coupling
tcoupl                   = v-rescale
nsttcouple               = 1
tc_grps                  = System
tau_t                    = 0.5
ref_t                    = 298.15
; Pressure coupling is on for NPT - we're doing NVT
pcoupl                   = no

; velocity generation
gen_vel                  = yes
gen-temp                 = 298.15
gen-seed                 = {randseed} ; need to randomize the seed each time.

; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
lincs-iter               = 2


; FREE ENERGY CONTROL OPTIONS =
free-energy   	        = expanded
calc-lambda-neighbors 	= -1
sc-alpha 		= 0.5    ;     0.5 
sc-power 		= 1      ;     keep this at 1 
sc-sigma 	        = 0.3    ;     0.5
couple-moltype 		= {couple_moltype}  ; ligand mol type
couple-lambda0 		= vdw-q
couple-lambda1 		= none
couple-intramol 	= yes
init-lambda-state	= {init_lambda_state}

nstexpanded             = 250   ; trying exchanges every 1 ps
nstdhdl                 = 250   ; dhdl snaps every 1 ps   <-- BINGO these must be set the same
dhdl-print-energy 	= total
nst-transition-matrix 	= 250000

lmc-seed                = {randseed} ; should be randomized
lmc-gibbsdelta          = -1 ; print all energies
symmetrized-transition-matrix = yes

lmc-stats                       = wang-landau
lmc-move                        = metropolized-gibbs
lmc-weights-equil               = wl-delta
weight-equil-wl-delta           = 0.00001
init-wl-delta                   = {wl_increment_in_kT}   ; in units kT -  MRS had 10.0 at first
separate-dhdl-file              = yes
wl-scale                        = 0.8
wl-ratio                        = 0.7

coul-lambdas         = {coul_lambdas_string}
vdw-lambdas          = {vdw_lambdas_string}
fep-lambdas          = {fep_lambdas_string}
init-lambda-weights  = {init_lambda_weights}

"""

        if not self.ligand_only:
            self.mdp_text +=  """; 
; pulling parameters
pull                     = umbrella
pull_ngroups             = 2
pull_ncoords             = 1
pull_group1_name         = {pull_group1_name}
pull_group2_name         = {pull_group2_name}
pull-geometry            = direction-periodic
pull_coord1_groups       = 1 2
pull-dim                 = Y Y Y
pull_coord1_rate         = 0.00
pull_coord1_k            = {pull_coord1_k}
pull-start               = no
pull-coord1-init         = {pull_coord1_init}
pull-nstxout             = 250   ; 1 ps
pull-nstfout             = 250   ; 1 ps
"""


        # format the mdpfile text!
        self.mdp_text = self.mdp_text.format(
            couple_moltype      = self.couple_moltype, 
            pull_group1_name    = self.pull_group1_name, 
            pull_group2_name    = self.pull_group2_name,
            pull_coord1_k       = self.pull_coord1_k,
            pull_coord1_init    = self.pull_coord1_init,
            fep_lambdas_string  = self.fep_lambdas_string,
            coul_lambdas_string = self.coul_lambdas_string,
            vdw_lambdas_string  = self.vdw_lambdas_string,
            init_lambda_state   = self.init_lambda_state,
            init_lambda_weights = self.init_lambda_weights_string,
            wl_increment_in_kT  = self.wl_increment_in_kT,
            randseed            = self.randseed )

        return self.header_desc + self.mdp_text

    def write_to_filename(self, filename):
        """Writes the ee mdpfile to specified filename."""

        fout = open(filename, 'w')
        fout.write( self.mdpfile_text() )         
        fout.close()


