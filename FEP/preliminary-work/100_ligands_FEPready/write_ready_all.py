#! /usr/bin/env python
for run in range(1,99):
    print("python make_fep_ready.py ../100_ligands/RUN%d ./RUN%d"%(run,run) )
    print("# cleanup!")
    print("rm ./#* ;  rm step*")
