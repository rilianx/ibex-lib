#!/bin/bash
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
cd $M
xargs -P 8 -L 1 $M/l3600coc.sh < $M/jobs_l3600coc.txt
xargs -P 8 -L 1 $M/abls.sh < $M/jobs_abls.txt
echo "cadena_sis: l3600coc $(ls $M/l3600coc | wc -l) abls $(ls $M/abls | wc -l)"
