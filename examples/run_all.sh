#{hc4,acidhc4,3bcidhc4,3bcidobjhc4,mohc,acidmohc,3bcidmohc}
#{art,compo,xn,xn_onlyy}

./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.00 > outputs/20141016-step0.00
./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.00 >> outputs/20141016-step0.00

#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.05 > outputs/20141016-step0.05
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.05 >> outputs/20141016-step0.05
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.1 > outputs/20141016-step0.10
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.1 >> outputs/20141016-step0.10
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.0 0.0 > outputs/20141016-reference
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.0 0.0 >> outputs/20141016-reference
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.15 > outputs/20141016-step0.15
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.15 >> outputs/20141016-step0.15
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.2 > outputs/20141016-step0.20
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.2 >> outputs/20141016-step0.20
