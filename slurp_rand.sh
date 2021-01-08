#!/bin/bash

#SBATCH -t 0-00:10:30
#SBATCH -N 3
#SBATCH --ntasks-per-node=10

cd /home/bisselin/Students_20/AnilMichael/

for n in 1000
do
  bsprun -n 1 ./sequential_output "$n" >> seq_times.json
  for N in 1 2 3 4 5 6 7 8 9
  do
    bsprun -n $((N*N)) ./parallel_output "$N" "$N" "$n" >> par_times.json
  done
done
