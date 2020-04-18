#!/usr/bin/env python

from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from simtk import openmm, unit
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
import parmed
import sys, subprocess

water_model = 'tip3p'
solvent_padding = 10.0 * unit.angstrom
ionic_strength = 100 * unit.millimolar # 100
pressure = 1.0 * unit.atmospheres
collision_rate = 91.0 / unit.picoseconds
temperature = 300.0 * unit.kelvin
timestep = 2.0 * unit.femtoseconds
xtc_step = 1000 # test
nsteps_equil = 500000 # test

protein_forcefield = 'amber14/protein.ff14SB.xml'

# trying to use SMIRNOFF causes a torsion error during SystemGenerator call
#small_molecule_forcefield = 'openff_unconstrained-1.1.0.offxml'
small_molecule_forcefield = 'openff-1.1.0'
#small_molecule_forcefield = 'gaff-2.11'
solvation_forcefield = 'amber14/tip3p.xml'

receptor_file = './fixed_receptor.pdb'

for ligand_ndx in range(int(sys.argv[1]),int(sys.argv[2])): # range(1,101):
    print(f'Processing RUN{ligand_ndx}')
    # create the system, solvate, and write integrator
    output_prefix = f'RUN{ligand_ndx}'
    ligand_file = f'{output_prefix}/LIG{ligand_ndx}_h.sdf'

    ligand = Molecule.from_file(ligand_file)
    ligand.to_file(f'{output_prefix}/LIG{ligand_ndx}_h.pdb', file_format='pdb')

    receptor = app.PDBFile(receptor_file)
    receptor_structure = parmed.load_file(receptor_file)
    ligand_structure = parmed.load_file(f'{output_prefix}/LIG{ligand_ndx}_h.pdb')
    complex_structure = receptor_structure + ligand_structure

    barostat = openmm.MonteCarloBarostat(pressure, temperature)
    forcefield_kwargs = {'removeCMMotion': False, 'ewaldErrorTolerance': 5e-04,
        'nonbondedMethod': app.PME, 'constraints': None, 'rigidWater': False}
    system_generator = SystemGenerator(forcefields=[protein_forcefield,solvation_forcefield],
        barostat=barostat, forcefield_kwargs=forcefield_kwargs, molecules=[ligand],
        small_molecule_forcefield=small_molecule_forcefield)
    
    modeller = app.Modeller(complex_structure.topology, complex_structure.positions)
    modeller.addSolvent(system_generator.forcefield, model='tip3p',
        padding=solvent_padding, ionicStrength=ionic_strength)
    
    system = system_generator.create_system(modeller.topology)
    solvated_structure = parmed.openmm.load_topology(modeller.topology,
        system, xyz=modeller.positions)

    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    with open(f'{output_prefix}/integrator.xml', 'w') as f:
        f.write(openmm.XmlSerializer.serialize(integrator))


    # minimize and equilibrate
    print('Minimizing...')
    platform = openmm.Platform.getPlatformByName('CUDA')
    platform.setPropertyDefaultValue('CudaDeviceIndex', sys.argv[3])
    context = openmm.Context(system, integrator, platform)
    context.setPositions(modeller.positions)
    openmm.LocalEnergyMinimizer.minimize(context)
    
    print('Equilibrating...')
    integrator.step(nsteps_equil)
    
    # save equilibrated pdb/state/system
#    with open(f'{output_prefix}/equilibrated_complex.pdb', 'w') as f:
#        app.PDBFile.writeFile(modeller.topology, context.getState(
#            getPositions=True,enforcePeriodicBox=True).getPositions(),
#            file=f, keepIds=True)

    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
#    with open(f'{output_prefix}/state.xml','w') as f:
#        f.write(openmm.XmlSerializer.serialize(state))

    system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
#    with open(f'{output_prefix}/system.xml','w') as f:
#        f.write(openmm.XmlSerializer.serialize(system))

    parmed_system = parmed.openmm.load_topology(modeller.topology,
        system, xyz=modeller.positions)
    parmed_system.save(f'{output_prefix}/conf.gro', overwrite=True)
    parmed_system.save(f'{output_prefix}/topol.top', overwrite=True)

## GMX CODE
#    try:
#        min_grompp_cmd = f'gmx grompp -f mdp/ee_minim.mdp -c {output_prefix}/conf.gro -p {output_prefix}/topol.top -o {output_prefix}/em.tpr -maxwarn 1'
#        print(min_grompp_cmd)
#        min_mdrun_cmd = f'cd {output_prefix}; gmx mdrun -v -deffnm em -gpu_id {sys.argv[3]}; cd ..'
#        print(min_mdrun_cmd)
#        equil_grompp_cmd = f'gmx grompp -f mdp/npt.mdp -c {output_prefix}/em.gro -p {output_prefix}/topol.top -o {output_prefix}/equil.tpr -maxwarn 1'
#        print(equil_grompp_cmd)
#        equil_mdrun_cmd = f'cd {output_prefix}; gmx mdrun -v -deffnm equil -gpu_id {sys.argv[3]}; cd ..'
#        print(equil_mdrun_cmd)
#        subprocess.check_output(min_grompp_cmd, shell=True)
#        subprocess.check_output(min_mdrun_cmd, shell=True)
#       subprocess.check_output(equil_grompp_cmd, shell=True)
#        subprocess.check_output(equil_mdrun_cmd, shell=True)
#    except Exception as e:
#        print(e)
#        continue
