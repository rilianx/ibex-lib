#{hc4,acidhc4,3bcidhc4,3bcidobjhc4,mohc,acidmohc,3bcidmohc}
#{art,compo,xn,xn_onlyy}

#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.10 > outputs/20141104-step0.1
./run_ncxp2.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.10 >> outputs/20141104-step0.1

./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.05 > outputs/20141104-step0.05

#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.00 > outputs/20141104-step0.00
./run_ncxp2.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.00 >> outputs/20141104-step0.00

#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.01 > outputs/20141104-step0.01
./run_ncxp2.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.01 >> outputs/20141104-step0.01

#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.15 > outputs/20141104-step0.15
./run_ncxp2.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.15 >> outputs/20141104-step0.15

#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.1 0.0 > outputs/20141104-alpha0.1
./run_ncxp2.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.1 0.0 >> outputs/20141104-alpha0.1

./run_ncxp2.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.05 >> outputs/20141104-step0.05

#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.3 0.0 > outputs/20141104-alpha0.3
#~ ./run_ncxp2.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.3 0.0 >> outputs/20141104-alpha0.3
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.7 0.0 > outputs/20141104-alpha0.7
#~ ./run_ncxp2.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.7 0.0 >> outputs/20141104-alpha0.7
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.9 0.0 > outputs/20141104-alpha0.9
#~ ./run_ncxp2.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.9 0.0 >> outputs/20141104-alpha0.9

#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.05 > outputs/20141104-step0.05
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.05 >> outputs/20141104-step0.05
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.1 > outputs/20141104-step0.10
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.1 >> outputs/20141104-step0.10
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.0 0.0 > outputs/20141104-reference
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.0 0.0 >> outputs/20141104-reference
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.15 > outputs/20141104-step0.15
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.15 >> outputs/20141104-step0.15
#~ 
#~ ./run_ncxp.sh "../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05 optim 3600 no acidhc4 0.8 compo 3 8 0.5 0.2 > outputs/20141104-step0.20
#~ ./run_ncxp.sh "../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05 optim 3600 ampl acidhc4 0.8 compo 3 8 0.5 0.2 >> outputs/20141104-step0.20
