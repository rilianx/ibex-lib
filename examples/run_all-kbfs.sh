#{hc4,acidhc4,3bcidhc4,3bcidobjhc4,mohc,acidmohc,3bcidmohc}
#{art,compo,xn,xn_onlyy}


./run_ncxp.sh "../../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05-mindepth optim 3600 no acidhc4 0.8 compo 3 -2 8 > outputs/20141027-LB_ubtree_mindepth-BS2
./run_ncxp.sh "../../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05-mindepth optim 3600 ampl acidhc4 0.8 compo 3 -2 8 >> outputs/20141027-LB_ubtree_mindepth-BS2

./run_ncxp.sh "../../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05-mindepth optim 3600 no acidhc4 0.8 compo 3 -5 8 > outputs/20141027-LBmindepth-BS5
./run_ncxp.sh "../../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05-mindepth optim 3600 ampl acidhc4 0.8 compo 3 -5 8 >> outputs/20141027-LBmindepth-BS5

./run_ncxp.sh "../../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05-mindepth-noubtree optim 3600 no acidhc4 0.8 compo 3 -1 8 > outputs/20141027-LB_mindepth-Gr
./run_ncxp.sh "../../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05-mindepth-noubtree optim 3600 ampl acidhc4 0.8 compo 3 -1 8 >> outputs/20141027-LB_mindepth-Gr

./run_ncxp.sh "../../../benchmarks/coconut-library1-noTOacidhc4/*.bc*" solveroptim-araya05-maxdepth-noubtree optim 3600 no acidhc4 0.8 compo 3 -1 8 > outputs/20141027-LB_maxdepth-Gr
./run_ncxp.sh "../../../benchmarks/coconut-library2-noTOacidhc4/*.nl" solveroptim-araya05-maxdepth-noubtree optim 3600 ampl acidhc4 0.8 compo 3 -1 8 >> outputs/20141027-LB_maxdepth-Gr

