#!/bin/bash
# job.sh <kind> <bch>: banco de pasadas, lado caja = SoPlex 4.0.2 (R2) o SoPlex 8.1 (con y sin escalado)
H=/home/iaraya/ibex2024/results/escalera_2026-11-06; S=$H/work/S100; M=$H/work/M
D=$(cd "$(dirname "$0")" && pwd)
cd /home/iaraya/ibex2024/examples/dfb; ulimit -v 6000000
k=$1; n=$(basename $2 .bch)
case $k in
  cmp) c=$S/bancocmp/sol.$n.txt; pr=sol; x="DFBH_REPLAY_LIN=compo"; tol=1e-6;;
  opt) c=$S/bancoopt/opt.$n.txt; [ -s $c ] || c=$M/banco/opt.$n.txt; pr=opt; x=""; tol=1e-9;;
  sol) c=$M/banco/sol.$n.txt; pr=sol; x=""; tol=1e-6;;
esac
[ -s $c ] || exit 0
SC="DFB_SX_BFRT=0 DFB_SX_DEVEX=0 DFBH_LPCAP=100"
dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol DFBH_CAJA_TOPE=100"
for v in s4 s8 s8n; do
  o=$D/out/$k.$n.$v.txt; [ -s $o ] && continue
  case $v in s4) e="";; s8) e="DFBH_CAJA_SOPLEX7=1 CAJA7_SIMP=0";; s8n) e="DFBH_CAJA_SOPLEX7=1 CAJA7_SIMP=0 CAJA7_SCALER=0";; esac
  env $x $SC $dir $e DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 ./ibexreplay_s8 $2 $pr $c 2>&1 | grep -E 'abcaja|replay|Assert|terminate' > $o
done
