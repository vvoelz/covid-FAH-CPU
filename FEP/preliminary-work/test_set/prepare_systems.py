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



for i in [4]: #range(84,101):
    path = f'LIG{i}'
    print(f'Processing LIG{i}')

#    try:
    if 1:
        # DO LIGAND THINGS
        try:
            ligand_off_molecule = Molecule(f'{path}/results{i}.sdf')
        except Exception as e:
            print(e)
            continue
        ligand_pdbfile = PDBFile(f'{path}/results{i}.pdb')
        force_field = ForceField('openff_unconstrained-1.1.0.offxml') #smirnoff99Frosst.offxml') #'openff-1.0.0.offxml')
        ligand_system = force_field.create_openmm_system(ligand_off_molecule.to_topology())
        ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology,
                                                        ligand_system,
                                                        xyz=ligand_pdbfile.positions)
    if 1:
        # DO PROTEIN THINGS
        receptor_file = 'receptor.pdb'
        fixed_receptor_file = f'{path}/fixed_receptor.pdb'
        omm_forcefield = app.ForceField('amber14-all.xml')
        fixer = PDBFixer(receptor_file) #filename='receptor.pdb')
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
        complex_structure.save(f'{path}/complex.pdb', overwrite=True)
        with open(f'{path}/system.xml', 'w') as f:
            f.write(XmlSerializer.serialize(complex_system))
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
        complex_structure.save(f'{path}/conf.gro', overwrite=True)
        complex_structure.save(f'{path}/topol.top', overwrite=True)
        print('Equilibrating...')
        simulation.step(10000)
        simulation.reporters.append(app.DCDReporter(f'{path}/traj{i}.dcd', 10000))
        simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, 
            potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
            speed=True, totalSteps=1000, separator='\t'))

        print('Running Production...')
        simulation.step(100000)
        print('Done!')

#    except Exception as e:
#        print(e)
#        pass


# In[ ]:





# In[ ]:





# In[ ]:




