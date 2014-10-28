#{hc4,acidhc4,3bcidhc4,3bcidobjhc4,mohc,acidmohc,3bcidmohc}
#{art,compo,xn,xn_onlyy}

./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-lbub optim 3600 no acidhc4 0.8 compo 3 -1 8 > outputs/20141028-LB_mindepth-lbub
./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-lbub optim 3600 ampl acidhc4 0.8 compo 3 -1 8 >> outputs/20141028-LB_mindepth-lbub

./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-gBFS optim 3600 no acidhc4 0.8 compo 3 -1 8 > outputs/20141028-LB_mindepth-GBFS
./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-gBFS optim 3600 ampl acidhc4 0.8 compo 3 -1 8 >> outputs/20141028-LB_mindepth-GBFS


