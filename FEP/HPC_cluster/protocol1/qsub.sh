#!/bin/sh
#PBS -l walltime=24:00:00   # max is 48 hrs on owlsnest2
#PBS -N GPU1
#PBS -q gpu        # change this to "large" if used on owlsnest2
#PBS -l nodes=1:ppn=2  # max is 4 on cb2rr but 2 on owlsnest2
#PBS -o GPU1 
#PBS 

cd $PBS_O_WORKDIR

. ~/.bashrc       # source my own bash so I can get rid of "conda init" error
conda activate openff     # my openmm specific conda environment
module load cuda/10.0.130   # load cuda/10 if used on owlsnest2


export CUDA_VISIBLE_DEVICES=0
python simulate.py  &

export CUDA_VISIBLE_DEVICES=1
python simulate2.py  &

# can list more if you want, but max is 4 on cb2rr and 2 on owlsnest2

wait

