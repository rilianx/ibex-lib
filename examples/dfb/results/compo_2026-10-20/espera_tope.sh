#!/bin/bash
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
while ps -eo args | grep -q '[s]onda_filas\.sh /'; do sleep 15; done
cd $M && xargs -P 2 -L 1 $M/sonda_tope.sh < $M/jobs_sonda.txt
echo "tope: $(ls $M/sonda_compo | grep -c _tope)"
