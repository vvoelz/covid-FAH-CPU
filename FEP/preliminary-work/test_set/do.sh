#!/bin/bash

for i in {1..100}; do
    mkdir LIG$i
    cp results$i.mol2 LIG$i
    cp results$i.pdb LIG$i
    obabel -imol2 LIG$i/results$i.mol2 -osdf > LIG$i/results$i.sdf; done
