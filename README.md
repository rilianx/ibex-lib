[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

ibex-lib
========

http://www.ibex-lib.org


Installation
------------
./waf configure --with-optim --with-ampl --with-affine --prefix=. --gaol-dir= --lp-lib=soplex

./waf install


Test Example
------------
__build__/plugins/optim/ibexopt plugins/optim/benchs/medium/ex7_3_4.bch --trace --lmode=abst


