ParamEstim
======
This is a README for knowledge on how to use this plugin on ibex solver.


## Configuration
### Default build
Must change "source" variable to set example to execute on file wscript on folder "param-estim".

`source = bld.path.ant_glob ("examples/**/robust_estim1.cpp")`

Choices for source are:

* robust_estim1.cpp

* robust_estim2.cpp

* robust_estim3.cpp

* robust_estim4.cpp

* estim_optim.cpp

### ibex configuration
Must add the param-estim flag to use it and as optional the default optimizer (optim).

`./waf configure --with-optim --with-param-estim --prefix=. --gaol-dir= --lp-lib=soplex --with-affine`

## Compile

`./waf install` cross fingers and hope for the best.

## Launch

### Arguments
* filename (string): house40
* epsilon (double): 1.e-3
* Q (int): 21
* prebc (double)
* precd (double)
* flist (int)
* gaplimit (int)
* nbrand (int)
* bisect method (string): rr or rr2
* random seed

### Execute

`__build__/plugins/param-estim/estim_optim house40 1.e-3 21 1.e-4 1.e-4 1 0 0 rr2 1`


### Extras
For more examples, you must copy the input file from the examples folder to `__build__`

`cp examples/estimparam/house40 __build__/plugins/param-estim/house40`

Replace "house40" with the filename you want to test.

