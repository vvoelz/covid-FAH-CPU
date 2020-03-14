# FEP scripts for FAH

In this directory are scripts for creating and building new GRO_A7 (gmx 5.0.4) tprs after each returned WU.

* `create_ee_mdp.py` - a script to create a new FEP project 
* `continue_a_tpr.py` - a scripta analysis

Also:
* `expanded.py`  - contains the  `expanded_ensemble_mdpfile()` class 

**IMPORTANT!!!** DO NOT mess with the organization of these scripts in the repo!!! 
This git repo is pulled onto the FAH servers, and the FAH work serve calles them directly!!!!!

## How to set up a project

The FEP project _MUST_ be prepared with standard naming conventions!!
(You can see some examples of this in [../ee-protein-ligand](../ee-protein-ligand) and [../ee-ligand-only](../ee-ligand-only) )

* the Gromacs topology `*.top` MUST have a moleculetype `LIG`
* the index file `*.ndx` MUST have an index group `LIG` 
* it's a very good idea to have the ligand residue be named `LIG` in the `*.gro` file

### Protein-ligand simulations ###

For a protein-ligand simulation, there is a tether between the ligand to the protein,
and we need to have prepared an index file with the atom groups `a1-Protein`, `a2-Ligand` and `Restraint-Distance`,
define as in the following example:
```
            ....
            [ a1-Protein ]
            678
            [ a2-Ligand ]
            1564
            [ Restraint-Distance ]
            678 1564
```
Moreover, there needs to be atom groups `Protein` and `non-Protein`, which the thermostat will control separately.

### Ligand-only simulations ###
        There needs to be atom groups `Water` and `non-Water`, which the thermostat will control separately.
