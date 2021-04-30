#!/bin/bash

#PBS -l select=1:ncpus=1:mpiprocs=1,place=free:exclhost
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR

mpiicc -O3 -std=gnu99 *.c

mpirun -machinefile $PBS_NODEFILE -np 1 ./a.out

#mpirun -trace -machinefile $PBS_NODEFILE -np 16 ./a.out
