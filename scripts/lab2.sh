#!/bin/bash

#PBS -l select=1:ncpus=12:mem=400m:ompthreads=12,place=free:exclhost
#PBS -l walltime=00:15:00
#PBS -q S3077545

cd $PBS_O_WORKDIR
#echo "I run on node: `uname -n`"
echo "My working directory is: $PBS_O_WORKDIR"
#echo "Assigned to me nodes are:"
cat $PBS_NODEFILE

icc -openmp -O3 -o prog_naive -std=gnu99  main_naive.c matrix_omp.c data.c
icc -openmp -O3 -o prog_2 -std=gnu99  main_parallel.c matrix_parallel.c data.c
for i in {1..12}
do
		echo "Threads: $i"
			echo "1st variant"
				OMP_NUM_THREADS=$i ./prog_naive
					echo "2nd variant"
						OMP_NUM_THREADS=$i OMP_SCHEDULE=static ./prog_2
					done

					num_threads=4
					for i in {static,auto,dynamic,dynamic\,8,guided}
					do
							echo "Threads: $num_threads"
								echo "Schedule: ($i)"
									OMP_NUM_THREADS=$num_threads OMP_SCHEDULE=$i ./prog_2
								done
