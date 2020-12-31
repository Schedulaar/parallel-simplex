#!/bin/bash

bspcxx sequential-simplex.cpp -o seq_output --std=c++11 -O3
bspcxx parallel-simplex.cpp -o par_output --std=c++11 -O3


for n in 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000
do
  ./seq_output "$n" >> seq_times.json
  for N in 1 2 3 4
  do
    ./par_output "$N" "$N" "$n" >> par_times.json
  done

done
