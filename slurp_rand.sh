#!/bin/bash

#SBATCH -t 0-00:10:30
#SBATCH -N 3
#SBATCH --ntasks-per-node=24

module load 2019
module load BSPonMPI/1.1-foss-2018b

cd /home/bisselin/Students20/AnilMichael/parallel-simplex

for n in 1000
do
  bsprun -n 1 ./sequential_output rand "$n" >> seq_times.json
  for N in 1 2 3 4 5 6 7 8 9
  do
    bsprun -n $((N*N)) ./parallel_output "$N" "$N" rand "$n" >> par_times.json
  done
done
