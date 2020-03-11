#!/usr/bin/env python
# coding: utf-8

# In[41]:


#get_ipython().run_line_magic('ls', '5a2i/')


# In[44]:


import parmed
from pdbfixer import PDBFixer
from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from simtk import unit
from simtk import openmm as mm
from simtk.openmm import app, XmlSerializer
from simtk.openmm.app import NoCutoff, HBonds, PDBFile
from sys import stdout
import os,glob

import warnings
warnings.filterwarnings("ignore")


# In[ ]:



for i in range(1,8):
    path = f'HIT{i}'
    print(f'Processing HIT{i}')

    try:
        # DO LIGAND THINGS
        ligand_off_molecule = Molecule(f'{path}/hit{i}.sdf')
        ligand_pdbfile = PDBFile(f'{path}/hit{i}.pdb')
        force_field = ForceField('openff-1.0.0.offxml')
        ligand_system = force_field.create_openmm_system(ligand_off_molecule.to_topology())
        ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology,
                                                        ligand_system,
                                                        xyz=ligand_pdbfile.positions)
        # DO PROTEIN THINGS
        receptor_file = f'{path}/receptor.mol2'
        fixed_receptor_file = f'{path}/fixed_receptor.pdb'
        omm_forcefield = app.ForceField('amber14-all.xml')
        fixer = PDBFixer(filename=receptor_file)        
        missingresidues = fixer.findMissingResidues()
        rezez = fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        missingatoms = fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_receptor_file, 'w'))
        fixed_receptor = PDBFile(fixed_receptor_file)
        receptor_system = omm_forcefield.createSystem(fixed_receptor.topology)
        receptor_structure = parmed.openmm.load_topology(fixed_receptor.topology,
                                                         receptor_system,
                                                         xyz=fixed_receptor.positions)
        complex_structure = receptor_structure + ligand_structure
        complex_system = complex_structure.createSystem(nonbondedMethod=NoCutoff,
                                                        nonbondedCutoff=9.0*unit.angstrom,
                                                        constraints=HBonds,
                                                        removeCMMotion=False)
        complex_structure.save(f'{path}/hit{hit}_complex.pdb', overwrite=True)
        with open(f'{path}/hit{hit}_complex.xml', 'w') as f:
            f.write(XmlSerializer.serialize(complex_system))
        complex_structure.save(f'{path}/hit{hit}_complex.gro', overwrite=True)
        complex_structure.save(f'{path}/hit{hit}_complex.top', overwrite=True)
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 
            2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)

        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        simulation = app.Simulation(complex_structure.topology, complex_system, integrator, platform, 
            properties)
        simulation.context.setPositions(complex_structure.positions)

        print('Minimizing...')
        simulation.minimizeEnergy()

        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        print('Equilibrating...')
        simulation.step(10000)
        simulation.reporters.append(app.DCDReporter(f'{path}/hit{hit}_trajectory.dcd', 1000))
        simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, 
            potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
            speed=True, totalSteps=1000, separator='\t'))

        print('Running Production...')
        simulation.step(100000)
        print('Done!')

    except Exception as e:
        print(e)
        pass


# In[ ]:





# In[ ]:





# In[ ]:




