MAXSATZ
======
This plugin was created based on Maxsatz: A branch and bound solver for Max-SAT found in https://home.mis.u-picardie.fr/~cli/EnglishPage.html
Thanks to Chu Min LI, professor in University of Picardie Jules Verne.

This is a README knowledge on how to use this plugin on ibex solver.


## Configuration

### Default build

### ibex configuration

Must add the maxsat flag to use it and as optional the default optimizer (optim).

`./waf configure --with-optim --with-maxsat --prefix=. --gaol-dir= --lp-lib=soplex --with-affine`

## Compile

`./waf install` cross fingers and hope for the best.

## Launch

### Arguments

* filename (string): test_name (folder plugins/maxsat/test/)

### Execute

`__build__/plugins/maxsat/estim_optim plugins/maxsat/test/test_name`

Replace "test_name" with the filename you want to test.