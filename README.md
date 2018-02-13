[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

ibex-lib
========

http://www.ibex-lib.org

Instalation
-----------

Instalar g++, clang, bison, flex, zlib1g-dev

./waf configure --with-optim --with-pso --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex

./waf install

Para ejecutar:

 __build__/plugins/optim/ibexopt plugins/optim/benchs/coconutbenchmark-library2/dualc2.nl --trace
 
 __build__/plugins/pso/ibexpso plugins/optim/benchs/coconutbenchmark-library2/dualc2.nl --trace
 
Prueba sin restricciones
__build__/plugins/pso/ibexpso plugins/optim/benchs/benchs-unconstrainedoptim/dixon-price5.bch --trace

Prueba para configuracion de comportamiento
__build__/plugins/pso/ibexpso plugins/pso/benchs/cristian-price2.bch --trace

Graficar con matplotlib y python3
ejecutar pso con --trace para crear archivo con posiciones en raíz y luego ejecutar script
-> python3 plot.py