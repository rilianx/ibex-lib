#!/bin/bash
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
cd $M
xargs -P 8 -L 1 $M/d3600p.sh < $M/jobs_d3600p.txt
xargs -P 8 -L 1 $M/d7200p.sh < $M/jobs_d7200p.txt
xargs -P 8 -L 1 $M/prod.sh < $M/jobs_prod.txt
echo "cadena_larga: d3600p $(ls $M/d3600p | wc -l) d7200p $(ls $M/d7200p | wc -l) prod $(ls $M/prod | wc -l)"
