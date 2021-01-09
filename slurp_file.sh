#!/bin/bash

#SBATCH -t 0-00:10:30
#SBATCH -N 2
#SBATCH --ntasks-per-node=24

module load 2019
module load BSPonMPI/1.1-foss-2018b

cd /home/bisselin/Students20/AnilMichael/parallel-simplex

for file in 25FV47 ADLITTLE AFIRO AGG AGG2 AGG3 BANDM BEACONFD BLEND BNL1 BNL2 BRANDY D2Q06C DEGEN2 DEGEN3 E226 FFFFF800 ISRAEL LOTFI MAROS-R7 QAP8 QAP12 QAP15 SC50A SC50B SC105 SC205 SCAGR7 SCAGR25 SCFXM1 SCFXM2 SCFXM3 SCORPION SCRS8 SCSD1 SCSD6 SCSD8 SCTAP1 SCTAP2 SCTAP3 SHARE1B SHARE2B SHIP04L SHIP04S SHIP08L SHIP12L SHIP12S STOCFOR1 STOCFOR2 TRUSS WOOD1P WOODW
do
  bsprun -n 1 ./sequential_output file "netlib/$file" >> seq_times.json
  for N in 1 2 3 4 5 6
  do
    bsprun -n $((N*N)) ./parallel_output "$N" "$N" rand "netlib/$file" >> par_times.json
  done
done
