import os, sys
import time
import numpy as np

# Classes and Functions for preparing and continuing expanded-ensemble simulations

class expanded_ensemble_mdpfile(object):

    def __init__(self,  couple_moltype = '1MQ',
			pull_group1_name = 'a1-Protein', 
                        pull_group2_name = 'a2-Ligand', 
                        pull_coord1_k    = 200.0, 
                        pull_coord1_init  = 0.4, 
                        fep_lambdas      = np.arange(0.00, 1.05, 0.05),
                        init_lambda_weights = np.zeros(21),
                        init_lambda_state   = 20,
                        wl_increment_in_kT  = 3.0 ):

        """Initialize the class."""

        self.couple_moltype = couple_moltype 
        self.pull_group1_name = pull_group1_name
        self.pull_group2_name = pull_group2_name
        self.pull_coord1_k    = pull_coord1_k     # the spring constant in kJ/nm^2
        self.pull_coord1_init = pull_coord1_init  # the equilibrium distance between the pull group 1 and 2
        self.fep_lambdas      = fep_lambdas
        self.nlambdas         = len(self.fep_lambdas)
        self.fep_lambdas_string   = ' '.join(['%2.3f'%fep_lambdas[i] for i in range(self.nlambdas)])

        self.init_lambda_weights  = init_lambda_weights
        self.init_lambda_weights_string =  ' '.join(['%2.3f'%init_lambda_weights[i] for i in range(self.nlambdas)])

        self.init_lambda_state    = init_lambda_state
        self.wl_increment_in_kT   = wl_increment_in_kT


    def mdpfile_text(self):
        """Returns a string corresponding to the mdpfile contents."""

        self.header_desc = ";       generated from expanded_ensemble_mdpfile() on %s \n;\n;\n"%time.asctime()

        self.mdp_text = """; VARIOUS PREPROCESSING OPTIONS
; Preprocessor information: use cpp syntax.
; e.g.: -I/home/joe/doe -I/home/mary/roe
include                  = -I../top
; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
define                   = 

; RUN CONTROL PARAMETERS
integrator               = md-vv   ;  only md-vv works with ee;  sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 500000   ; 1000 ps = 1 ns
; For exact run continuation or redoing part of a run
init_step                = 0
; Part index is updated automatically on checkpointing (keeps files separate)
simulation_part          = 1
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 10
; group(s) for center of mass motion removal
comm_grps                = System

; LANGEVIN DYNAMICS OPTIONS
; Friction coefficient (amu/ps) and random seed
bd-fric                  = 0
ld-seed                  = -1 

; TEST PARTICLE INSERTION OPTIONS
rtpi                     = 0.05

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 500000  ; 1 ns 
nstvout                  = 500000  ; 1 ns
nstfout                  = 0
; Output frequency for energies to log file and energy file
nstlog                   = 500     ; 1 ps
nstcalcenergy            = 500     ; 1 ps
nstenergy                = 500     ; 1 ps
; Output frequency and precision for .xtc file
nstxtcout                = 500     ; every 1 ps
xtc-precision            = 1000
; This selects the subset of atoms for the .xtc file. You can
; select multiple groups. By default all atoms will be written.
xtc_grps                 = 1MQ
; Selection of energy groups
energygrps               = WaterProtein Water_and_ions 1MQ

; NEIGHBORSEARCHING PARAMETERS
; nblist update frequency
nstlist                  = 10
; ns algorithm (simple or grid)
ns_type                  = grid
; Periodic boundary conditions: xyz, no, xy
pbc                      = xyz
periodic_molecules       = no
; nblist cut-off        
rlist                    = 0.9
; long-range cut-off for switched potentials
rlistlong                = -1

; OPTIONS FOR ELECTROSTATICS AND VDW
; Method for doing electrostatics
coulombtype              = PME
rcoulomb-switch          = 0
rcoulomb                 = 0.9
; Relative dielectric constant for the medium and the reaction field
epsilon_r                = 1
epsilon_rf               = 1
; Method for doing Van der Waals
vdw-type                 = Cut-off
; cut-off lengths       
rvdw-switch              = 0
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = No
; Extension of the potential lookup tables beyond the cut-off
table-extension          = 1
; Seperate tables between energy group pairs
energygrp_table          = 
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; FFT grid size, when a value is 0 fourierspacing will be used
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters
pme_order                = 4
ewald_rtol               = 1e-5
ewald_geometry           = 3d
epsilon_surface          = 0
optimize_fft             = no

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling  
tcoupl                   = Berendsen
nsttcouple               = -1
nh-chain-length          = 10
; Groups to couple separately
tc-grps                  = Protein non-Protein
; Time constant (ps) and reference temperature (K)
tau_t                    = 1.0  1.0
ref_t                    = 300  300
; Pressure coupling     
Pcoupl                   = no
Pcoupltype               = Isotropic
nstpcouple               = -1
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                    = 1.0
compressibility          = 4.5e-5
ref_p                    = 1.0
; Scaling of reference coordinates, No, All or COM
refcoord_scaling         = No
; Random seed for Andersen thermostat
andersen_seed            = 815131

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                  = yes
gen_temp                 = 300
;gen_seed		 = -1 <<<--- FAHWorkServer complains that there are multiple lines; it adds its own line!

; OPTIONS FOR BONDS    
constraints              = hbonds
; Type of constraint algorithm
constraint-algorithm     = Lincs
; Do not constrain the start configuration
continuation             = no
; Use successive overrelaxation to reduce the number of shake iterations
Shake-SOR                = no
; Relative tolerance of shake
shake-tol                = 0.0001
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 4
; Number of iterations in the final step of LINCS. 1 is fine for
; normal simulations, but use 2 to conserve energy in NVE runs.
; For energy minimization with constraints it should be 4 to 8.
lincs-iter               = 1
; Lincs will write a warning to the stderr if in one step a bond
; rotates over more degrees than
lincs-warnangle          = 30
; Convert harmonic bonds to morse potentials
morse                    = no

; ENERGY GROUP EXCLUSIONS
; Pairs of energy groups for which all non-bonded interactions are excluded
energygrp_excl           = 

; WALLS                
; Number of walls, type, atom types, densities and box-z scale factor for Ewald
nwall                    = 0
wall_type                = 9-3
wall_r_linpot            = -1
wall_atomtype            = 
wall_density             = 
wall_ewald_zfac          = 3

; Free energy control stuff ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
free-energy              = expanded
init-lambda-state        = {init_lambda_state}
;  delta-lambda             = 0
fep-lambdas              = {fep_lambdas_string}
calc-lambda-neighbors    = -1
sc-alpha                 = 0
sc-power                 = 0
sc-sigma                 = 0.3
nstdhdl                  = 500
separate-dhdl-file       = yes
dhdl-derivatives         = yes
;   dh_hist_size             = 0
;   dh_hist_spacing          = 0.1
couple-moltype           = {couple_moltype}
couple-lambda0           = none ; off
couple-lambda1           = vdw-q ; on
couple-intramol          = no


; expanded ensemble variables
nstexpanded              = 500     ; 1 ps
lmc-stats                = wang-landau
lmc-move                 = metropolized-gibbs
lmc-weights-equil        = no
weight-equil-number-all-lambda = -1
weight-equil-number-samples = -1
weight-equil-number-steps = -1
weight-equil-wl-delta    = -1
weight-equil-count-ratio = -1

; Seed for Monte Carlo in lambda space
lmc-seed                 = -1
mc-temperature           = -1
lmc-repeats              = 1
lmc-gibbsdelta           = -1
lmc-forced-nstart        = 0
symmetrized-transition-matrix = no
nst-transition-matrix    = -1
mininum-var-min          = 100
init-lambda-weights      = {init_lambda_weights}
weight-c-range           = 0
wl-scale                 = 0.8
wl-ratio                 = 0.8
init-wl-delta            = {wl_increment_in_kT}  ; kT -- default is 1.0
wl-oneovert              = no

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
pull-nstxout             = 500   ; 1 ps
pull-nstfout             = 500   ; 1 ps
""".format(couple_moltype      = self.couple_moltype, 
           pull_group1_name    = self.pull_group1_name, 
           pull_group2_name    = self.pull_group2_name,
           pull_coord1_k       = self.pull_coord1_k,
           pull_coord1_init    = self.pull_coord1_init,
           fep_lambdas_string  = self.fep_lambdas_string,
           init_lambda_state   = self.init_lambda_state,
           init_lambda_weights = self.init_lambda_weights_string,
           wl_increment_in_kT  = self.wl_increment_in_kT )

        return self.header_desc + self.mdp_text

    def write_to_filename(self, filename):
        """Writes the ee mdpfile to specified filename."""

        fout = open(filename, 'w')
        fout.write( self.mdpfile_text() )         
        fout.close()



### A derived class for an mdpfile wth NO RESTRAINT
   
class expanded_ensemble_mdpfile_NOREST(expanded_ensemble_mdpfile):

    def __init__(self,  couple_moltype = '1MQ',
                        fep_lambdas      = np.arange(0.00, 1.05, 0.05),
                        init_lambda_weights = np.zeros(21),
                        init_lambda_state   = 20,
                        wl_increment_in_kT  = 3.0 ):

        """Initialize the class."""

        self.couple_moltype = couple_moltype 
        self.fep_lambdas      = fep_lambdas
        self.nlambdas         = len(self.fep_lambdas)
        self.fep_lambdas_string   = ' '.join(['%2.3f'%fep_lambdas[i] for i in range(self.nlambdas)])

        self.init_lambda_weights  = init_lambda_weights
        self.init_lambda_weights_string =  ' '.join(['%2.3f'%init_lambda_weights[i] for i in range(self.nlambdas)])

        self.init_lambda_state    = init_lambda_state
        self.wl_increment_in_kT   = wl_increment_in_kT


    def mdpfile_text(self):
        """Returns a string corresponding to the mdpfile contents.

        NOTE -- these will be used for ligand-only simulations in small boxes of water, so total time is 10 ns"""

        self.header_desc = ";       generated from expanded_ensemble_mdpfile() on %s \n;\n;\n"%time.asctime()

        self.mdp_text = """; VARIOUS PREPROCESSING OPTIONS
; Preprocessor information: use cpp syntax.
; e.g.: -I/home/joe/doe -I/home/mary/roe
include                  = -I../top
; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
define                   = 

; RUN CONTROL PARAMETERS
integrator               = md-vv   ;  only md-vv works with ee;  sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 5000000   ; 10 000 ps = 10 ns
; For exact run continuation or redoing part of a run
init_step                = 0
; Part index is updated automatically on checkpointing (keeps files separate)
simulation_part          = 1
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 10
; group(s) for center of mass motion removal
comm_grps                = System

; LANGEVIN DYNAMICS OPTIONS
; Friction coefficient (amu/ps) and random seed
bd-fric                  = 0
ld-seed                  = -1 

; TEST PARTICLE INSERTION OPTIONS
rtpi                     = 0.05

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 500000  ; 1 ns 
nstvout                  = 500000  ; 1 ns
nstfout                  = 0
; Output frequency for energies to log file and energy file
nstlog                   = 5000     ; 10 ps
nstcalcenergy            = 5000     ; 10 ps
nstenergy                = 5000     ; 10 ps
; Output frequency and precision for .xtc file
nstxtcout                = 50000  ; every 100 ps
xtc-precision            = 1000
; This selects the subset of atoms for the .xtc file. You can
; select multiple groups. By default all atoms will be written.
xtc_grps                 = 1MQ
; Selection of energy groups
energygrps               = Water non-Water

; NEIGHBORSEARCHING PARAMETERS
; nblist update frequency
nstlist                  = 10
; ns algorithm (simple or grid)
ns_type                  = grid
; Periodic boundary conditions: xyz, no, xy
pbc                      = xyz
periodic_molecules       = no
; nblist cut-off        
rlist                    = 0.9
; long-range cut-off for switched potentials
rlistlong                = -1

; OPTIONS FOR ELECTROSTATICS AND VDW
; Method for doing electrostatics
coulombtype              = PME
rcoulomb-switch          = 0
rcoulomb                 = 0.9
; Relative dielectric constant for the medium and the reaction field
epsilon_r                = 1
epsilon_rf               = 1
; Method for doing Van der Waals
vdw-type                 = Cut-off
; cut-off lengths       
rvdw-switch              = 0
rvdw                     = 0.9
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = No
; Extension of the potential lookup tables beyond the cut-off
table-extension          = 1
; Seperate tables between energy group pairs
energygrp_table          = 
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; FFT grid size, when a value is 0 fourierspacing will be used
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters
pme_order                = 4
ewald_rtol               = 1e-5
ewald_geometry           = 3d
epsilon_surface          = 0
optimize_fft             = no

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling  
tcoupl                   = Berendsen
nsttcouple               = -1
nh-chain-length          = 10
; Groups to couple separately
tc-grps                  = Water   non-Water
; Time constant (ps) and reference temperature (K)
tau_t                    = 1.0  1.0
ref_t                    = 300  300
; Pressure coupling     
Pcoupl                   = no
Pcoupltype               = Isotropic
nstpcouple               = -1
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                    = 1.0
compressibility          = 4.5e-5
ref_p                    = 1.0
; Scaling of reference coordinates, No, All or COM
refcoord_scaling         = No
; Random seed for Andersen thermostat
andersen_seed            = 815131

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                  = yes
gen_temp                 = 300
;gen_seed		 = -1 <<<--- FAHWorkServer complains that there are multiple lines; it adds its own line!

; OPTIONS FOR BONDS    
constraints              = hbonds
; Type of constraint algorithm
constraint-algorithm     = Lincs
; Do not constrain the start configuration
continuation             = no
; Use successive overrelaxation to reduce the number of shake iterations
Shake-SOR                = no
; Relative tolerance of shake
shake-tol                = 0.0001
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 4
; Number of iterations in the final step of LINCS. 1 is fine for
; normal simulations, but use 2 to conserve energy in NVE runs.
; For energy minimization with constraints it should be 4 to 8.
lincs-iter               = 1
; Lincs will write a warning to the stderr if in one step a bond
; rotates over more degrees than
lincs-warnangle          = 30
; Convert harmonic bonds to morse potentials
morse                    = no

; ENERGY GROUP EXCLUSIONS
; Pairs of energy groups for which all non-bonded interactions are excluded
energygrp_excl           = 

; WALLS                
; Number of walls, type, atom types, densities and box-z scale factor for Ewald
nwall                    = 0
wall_type                = 9-3
wall_r_linpot            = -1
wall_atomtype            = 
wall_density             = 
wall_ewald_zfac          = 3

; Free energy control stuff ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
free-energy              = expanded
init-lambda-state        = {init_lambda_state}
;  delta-lambda             = 0
fep-lambdas              = {fep_lambdas_string}
calc-lambda-neighbors    = -1
sc-alpha                 = 0
sc-power                 = 0
sc-sigma                 = 0.3
nstdhdl                  = 5000
separate-dhdl-file       = yes
dhdl-derivatives         = yes
;   dh_hist_size             = 0
;   dh_hist_spacing          = 0.1
couple-moltype           = {couple_moltype}
couple-lambda0           = none ; off
couple-lambda1           = vdw-q ; on
couple-intramol          = no


; expanded ensemble variables
nstexpanded              = 5000     ; 10 ps
lmc-stats                = wang-landau
lmc-move                 = metropolized-gibbs
lmc-weights-equil        = no
weight-equil-number-all-lambda = -1
weight-equil-number-samples = -1
weight-equil-number-steps = -1
weight-equil-wl-delta    = -1
weight-equil-count-ratio = -1

; Seed for Monte Carlo in lambda space
lmc-seed                 = -1
mc-temperature           = -1
lmc-repeats              = 1
lmc-gibbsdelta           = -1
lmc-forced-nstart        = 0
symmetrized-transition-matrix = no
nst-transition-matrix    = -1
mininum-var-min          = 100
init-lambda-weights      = {init_lambda_weights}
weight-c-range           = 0
wl-scale                 = 0.8
wl-ratio                 = 0.8
init-wl-delta            = {wl_increment_in_kT}  ; kT -- default is 1.0
wl-oneovert              = no

; pulling parameters
pull                     = no
""".format(couple_moltype      = self.couple_moltype, 
           fep_lambdas_string  = self.fep_lambdas_string,
           init_lambda_state   = self.init_lambda_state,
           init_lambda_weights = self.init_lambda_weights_string,
           wl_increment_in_kT  = self.wl_increment_in_kT )

        return self.header_desc + self.mdp_text

    def write_to_filename(self, filename):
        """Writes the ee mdpfile to specified filename."""

        fout = open(filename, 'w')
        fout.write( self.mdpfile_text() )         
        fout.close()
   

if __name__ == "__main__": 
    e = expanded_ensemble_mdpfile()
    print(e.mdpfile_text())
    e.write_to_filename('test.mdp')

    e2 = expanded_ensemble_mdpfile_NOREST()
    print(e2.mdpfile_text())
    e2.write_to_filename('test-NOREST.mdp')

 



