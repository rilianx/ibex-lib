ibex-lib
========

http://www.ibex-lib.org

Instalation
-----------

Must install g++, clang, bison, flex, zlib1g-dev

On folder ibex-lib:
./waf configure --with-optim --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex

Then compile the solver using:
./waf install
     
To execute:
__build__/plugins/optim/ibexopt plugins/optim/benchs/coconutbenchmark-library2/dualc2.nl --trace