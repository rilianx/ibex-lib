ParamEstim
======
This is a README for knowledge on how to use this plugin on ibex solver.


## Configuration
### wscript
Must change "source" variable to set example to execute on file wscript on folder "param-estim".

`source = bld.path.ant_glob ("examples/**/robust_estim1.cpp")`

Choices for source (for now) are:

* robust_estim1.cpp

* robust_estim2.cpp

* robust_estim3.cpp

* robust_estim4.cpp

### ibex configuration
Must add the param-estim flag to use it and as optional the default optimizer (optim).

`./waf configure --with-optim --with-param-estim --prefix=. --gaol-dir= --lp-lib=soplex --with-affine`

### extras
For some examples, you must copy the input file from the examples folder to `__build__`

`cp examples/estimparam/house40 __build__/plugins/param-estim/house40`

Replace "house40" with the filename you want to test.

## Compile

`./waf install` cross fingers and hope for the best.

## Launch

`__build__/plugins/param-estim/estim_house_optim house40 1.e-3 21 1.e-4 1.e-4 1 0 0 rr2 1`