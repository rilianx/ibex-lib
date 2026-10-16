#!/bin/bash
# Banco unico de 60 s (2026-10-16): f60 (sistemas + 139 opt, wall con ibexopt_fix) + 101 opt nuevos
# (sc/caja/lejos desde h3600 y fix/h3600.*, el resto de brazos desde fix/f60.*), sin robot.
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
rm -rf $M/u60; mkdir $M/u60
cp $M/f60/*.csv $M/u60/
cp $M/h3600/*.csv $M/u60/
for f in $M/fix/h3600.*.csv; do b=${f##*/}; cp $f $M/u60/${b#h3600.}; done
for f in $M/fix/f60.*.csv; do b=${f##*/}; cp $f $M/u60/${b#f60.}; done
rm -f $M/u60/opt.robot.*
ls $M/u60 | grep '^opt' | cut -d. -f2 | sort -u | wc -l
ls $M/u60 | grep '^opt' | cut -d. -f3 | sort | uniq -c
