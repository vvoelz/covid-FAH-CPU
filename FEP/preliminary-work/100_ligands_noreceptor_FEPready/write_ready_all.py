#! /usr/bin/env python
for run in range(1,101):
    print("python ../scripts/make_fep_ready.py ../100_ligands_noreceptor/RUN%d ./RUN%d ligonly"%(run,run) )
    print("# cleanup!")
    print("rm ./#* ;  rm step*")
