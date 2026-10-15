#!/bin/bash
# Sonda de costo rehecha (2026-10-14): la caja negra con su propio orden (sistemas cerca, opt lejos), no el del tableau viejo de DFB.
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/costo2/$2.$n.txt
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $2 = sol ]; then
  DFBH_ORDEN_CAJA=cerca DFBH_ABCAJA=1 timeout 200 ./ibexsolve_dfb_oc "$1" --lr=dfbhull --lin=art -t 60 --seed=1 2>&1 | grep '\[abcaja\]' > $o
else
  DFBH_ORDEN_CAJA=lejos DFBH_ABCAJA=1 DFBH_CAJA_TOL=1e-9 timeout 250 ./ibexopt_oc "$1" --filtering=acidhc4 --lr=dfbhull --seed=1 -t60 2>&1 | grep '\[abcaja\]' > $o
fi
[ -s "$o" ] || echo NA > $o
