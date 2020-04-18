#!/usr/bin/env python

import itertools
import subprocess
import mdtraj as md

for i in [1]: #range(1,101):
    structure = md.load(f'RUN{i}/conf.gro')
    protein_anchor = [a.index for a in structure.topology.atoms if a.residue.index == 164 and a.name == 'CA']
    ligand_atoms = [a.index for a in structure.topology.atoms if a.residue.name == 'UNL' and a.element.symbol != 'H']
    pairs = [list(x) for x in itertools.product(protein_anchor, ligand_atoms)]
    distances = list(md.compute_distances(structure,pairs))
    ligand_anchor = ligand_atoms[distances.index(min(distances))]
    print(protein_anchor, ligand_anchor)
    #?????
    cmd = f'echo -e "a{protein_anchor[0] + 1}\na{ligand_anchor + 1}\nq\n" | gmx make_ndx -f RUN{i}/conf.gro -o RUN{i}/ndx.ndx'
    subprocess.check_output(cmd, shell=True)
    cmd = f'echo -e "name 24 lol" | gmx make_ndx -n RUN{i}/ndx.ndx -o RUN{i}/lol.ndx'
    subprocess.check_output(cmd, shell=True)
    #name 24 a1_Protein\na{ligand_anchor + 1}\nname 25 a2_Ligand\nq\n" | gmx make_ndx -f RUN{i}/conf.gro -o RUN{i}/index.ndx'
#    print(cmd)
#    subprocess.check_output(cmd, shell=True)
     
