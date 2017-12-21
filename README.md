[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

ibex-lib
========

http://www.ibex-lib.org

Instalation
-----------

Instalar g++, clang, bison, flex, zlib1g-dev

./waf configure --with-optim  --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex

./waf install

export PKG_CONFIG_PATH=/directorio_ibex/ibex-lib/share/pkgconfig   -> Esto se tiene que hacer cada vez que se conecta a la maquina

export PKG_CONFIG_PATH=/home/iaraya/github/ibex/ibex-dev-dag/ibex-lib/share/pkgconfig



Para compilar todo en una sola linea yo hago lo siguiente:
1. ingreso al directorio raiz de ibex (solo la primera vez)
2. cd plugins/optim/main (solo la primera vez)
3. cd -; sudo ./waf install; cd -; rm optimisolver; make optimsolver (cada vez que quiero re-compilar todo)

Para ejecutar:

 __build__/plugins/optim/ibexopt plugins/optim/benchs/coconutbenchmark-library2/dualc8.nl --trace
 
 __build__/plugins/pso/ibexpso plugins/optim/benchs/coconutbenchmark-library2/dualc8.nl --trace