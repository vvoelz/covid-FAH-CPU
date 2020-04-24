#!/bin/bash

for i in RUN*; do
    sed -i "s/nstxout                  = 5000    ; every 10 ps/nstxout                  = 50000    ; every 100 ps/" $i/production.mdp
    sed -i "s/nstvout                  = 5000/nstvout                  = 50000/" $i/production.mdp
    sed -i "s/nstlog                   = 500/nstlog                   = 5000/" $i/production.mdp
    sed -i "s/nstxtcout                = 500 ; every 1 ps/nstxtcout                = 5000 ; every 1 ps/" $i/production.mdp
    sed -i "s/nstenergy                = 500 ; every 1 ps/nstenergy                = 5000 ; every 1 ps/" $i/production.mdp; done
