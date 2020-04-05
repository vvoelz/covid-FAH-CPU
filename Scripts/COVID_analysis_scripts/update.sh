#!/bin/bash

# This script will scrape data from FAH projects using scrape.py which will spit out pkl files which can be fed into
# analysis.py to plot free energy as a function of ligand-receptor coupling. This will update pkl's if they already exist.
# Example inputs
#python scrape.py 14363 MS0323_RL_1-500 &
#python scrape.py 14364 MS0323_L_1-500 &

desc=72 # Description i.e. ligand st
strt=1  # Ligand Starting Number
end=500 # Ligand Ending Number
step=500 # Step size i.e. Runs/project

for i in {1..12}; do     # How many projectsone wants to update
        proj=14600       # starting project -1 if indexing starts at 1
	python scrape.py $((proj+i)) ${desc}_RL_${strt}_${end} &
	python scrape.py $((proj+i)) ${desc}_RL_${strt}_${end} &
	strt=$((strt + step))
	end=$((end + step));
done

