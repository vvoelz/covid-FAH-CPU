This directory shows a working example that more directly emulates the creation and continuation
of ee WUs using the "master" version of our FAH scripts that live here:

`~/git/covid-FAH-CPU/FEP/scripts`

An example of creating an ee `*.tpr` from  is in `./create`:

 ```
# First, we create a ligand-only ee *.mdp for the *.tpr we're about to build
python ../scripts/create_ee_mdp.py RUN0/npt.gro RUN0/topol.top RUN0/index.ndx RUN0/prod.mdp ligonly

# Then, we grompp and make the *.tpr for the WU -- the FAH server code will do this
mkdir job0
$GMXBIN/gmx grompp -c RUN0/npt.gro -f RUN0/prod.mdp -p RUN0/topol.top -n RUN0/index.ndx -o job0/frame0.tpr -po job0/mdout.mdp -maxwarn 1

echo "# To test the tpr, do this:"
echo "cd job0"
echo "nohup gmx mdrun -v -s frame0.tpr &" 
```


