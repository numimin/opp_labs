#!/bin/bash

#PBS -l select=2:ncpus=8:mem=400m:mpiprocs=8,place=free:exclhost
#PBS -l walltime=00:15:00

cd $PBS_O_WORKDIR
#echo "I run on node: `uname -n`"
echo "My working directory is: $PBS_O_WORKDIR"
#echo "Assigned to me nodes are:"
cat $PBS_NODEFILE

MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
#echo $MPI_NP

mpiicc -O3 -o prog_1 -std=gnu99  main.c matrix.c data.c split_by_row.c
mpiicc -O3 -o prog_2 -std=gnu99  main.c matrix.c data.c split_by_column.c
icc -O3 -o prog_naive -std=gnu99  main_naive.c matrix.c data.c
./prog_naive
mpirun -machinefile $PBS_NODEFILE -np 2 ./prog_2
mpirun -machinefile $PBS_NODEFILE -np 4 ./prog_2
mpirun -machinefile $PBS_NODEFILE -np 8 ./prog_2
mpirun -trace -machinefile $PBS_NODEFILE -np 16 ./prog_2
mpirun -machinefile $PBS_NODEFILE -np 2 ./prog_1
mpirun -machinefile $PBS_NODEFILE -np 4 ./prog_1
mpirun -machinefile $PBS_NODEFILE -np 8 ./prog_1
mpirun -trace -machinefile $PBS_NODEFILE -np 16 ./prog_1
