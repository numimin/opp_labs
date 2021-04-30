#!/bin/bash

#PBS -l select=2:ncpus=8:mem=2000m:mpiprocs=8,place=free:exclhost
#PBS -l walltime=00:15:00

cd $PBS_O_WORKDIR
echo "I run on node: `uname -n`"
echo "My working directory is: $PBS_O_WORKDIR"
echo "Assigned to me nodes are:"
cat $PBS_NODEFILE

mpiicc -O3 -o prog -std=gnu99  main.c data.c
for i in {1,2,4,6,8,12,16}
do
	mpirun -machinefile $PBS_NODEFILE -np $i ./prog
done

mpirun -trace -machinefile $PBS_NODEFILE -np 16 ./prog
