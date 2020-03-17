This folder includes all work done to get a protocol that people can use on HPC clusters (CB2RR/owlsnest2) to prepare the systems and run OpenMM simulations for minimization and equilibration. The output files are gro and top files that can be used for Gromacs production run. The scripts should be easily modified to have Gromacs minimization and equilibration done. Since running OpenMM simulations are more tricky than Gromacs, here I only focused on the OpenMM simulations. This work cannt be done without Matt's help.  --Yunhui 03/2020

First let's set up your CB2RR/owlsnest2 conda environment.

1. Check your condarc file:

```cat ~/.condarc```

You should see these channels:
```channels:
  - conda-forge
  - defaults
  - omnia
  - omnia-dev
```

2. create a new environment for OpenMM 7.5

```conda create -n openff python==3.7 openmm==7.5.0 openforcefield```

```conda activate openff```

3. Normally step 2 will install OpenMM 7.5 with CUDA 10.2 build. However, on CB2RR/owlsnest2, only CUDA 10.0 is available. 
So let's re-install OpenMM with CUDA 10.0 build. 

```conda install -c omnia-dev/label/cuda100 openmm==7.5.0``` (Thanks John's help here)

4. We also need to install ```openmmforcefields```:

```conda install --yes -c conda-forge -c omnia openmmforcefields```

Okay, now you should have everything ready on CB2RR/owlsnest2.

Next, let's move on to the running jobs part.

There are two protocols available. Both are doing the same thing but in a different way. I think they both might be useful in different cases so I have both uploaded here.

Both protocol requires folders like "RUN1, RUN2, etc" for each ligand and in each folder there should be a pdb file and two sdf files. One sdf file was originally given by our collaborators and Dylan's scripts will do protonation work and generate another sdf file and a pdb file. (The scripts can be found here:https://github.com/vvoelz/covid-FAH-CPU/blob/master/Scripts/COVID-19_SDF_parser.py)
Also a pdb file of the receptor is required -- "fixed_receptor.pdb".

Protocol1:
The python script "simulate.py" has a "for loop" for multiple ligands.
People can directly run simulate.py with desired number of processors (check qsub.sh for details) by preparing multiple "simulate.py"-like files for each processor (GPU). The only thing need to be highlighted is that in the end of the "simulate.py" you can see several "del" lines. This is very important and necessary to run this script on HPC clusters (CB2RR tested, owlsnest2 not tested yet).

Protocol2:
The python script "simulate.py" is designed for only 1 ligand with a sys.argv for the ligand index. The runme file will have a "for loop" for mutiple ligands. The difference is in the "simulate.py", you don't need "del" lines in the end. Correspondingly, in the qsub file, you are running "runme" file instead of "python" files used in protocol1.

Both protocol should have similar performance. My test is still running and I will update once they are done. The current estimate is ~ 20 min for 1 ligand (prep+min+equil). That said, on CB2RR we should have ~ 280 ligands finished by 1 GPU node per day.
