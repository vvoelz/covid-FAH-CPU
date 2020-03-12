#!/bin/bash

for i in $(find . | grep topol | sed "s/.topol.*//"); do
    sed -i "s/\[ atomtypes \]/\n#include \"amber99sb-ildn.ff\/forcefield.itp\"\n\n\[ atomtypes \]/" $i/topol.top
    sed -i  "s/\[ molecules \]/\n#include \"amber99sb-ildn.ff\/tip3p.itp\"\n\n#include \"amber99sb-ildn.ff\/ions.itp\"\n\n\[ molecules \]/" $i/topol.top
    cd $i
    gmx editconf -f conf.gro -bt cubic -d 1.5 -o box.gro
    gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro &&
    gmx grompp -f ../ions.mdp -c solv.gro -p topol.top -o ions.tpr &&
    gmx genion -s ions.tpr -pname NA -nname CL -neutral -conc 0.1 -o ions.gro &&
    gmx grompp -f ../minim.mdp -c ions.gro -p topol.top -o em.tpr &&
    gmx mdrun -v -deffnm em
    cd ..
    done
