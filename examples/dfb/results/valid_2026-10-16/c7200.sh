#!/bin/bash
# Caja a 7200 s en los 4 problemas que solo SC termina a 3600 s (hs101 hs112 ex6_2_13 pollut), para tener pares y no censura. $1 nombre $2 semilla
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
cd ~/ibex2024/examples/dfb
o=$M/c7200/opt.$1.caja.$2.csv
[ -s "$o" ] && exit 0
ulimit -v 7000000
r=$(env DFBH_ORDEN=lejos DFBH_MOTOR=caja DFBH_CAJA_TOL=1e-9 timeout 7400 ./ibexopt /home/iaraya/ibex2026/benchs/optim/all/$1.bch --filtering=acidhc4 --lr=dfbhull --seed=$2 -t7200 2>/dev/null | tail -1)
[ -z "$r" ] && r="NA"
echo "$1,caja,$2,$r" > $o
