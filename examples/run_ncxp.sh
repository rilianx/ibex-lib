
dir=$1
algo=$2
type=$3
time=$4
ampl=$5


filt=$6  #{hc4,acidhc4,3bcidhc4,3bcidobjhc4,mohc,acidmohc,3bcidmohc}
taumohc=$7
lr=$8  #{art,compo,xn,xn_onlyy}
rand=$9
processors=${10}
alpha=${11}
step=${12}


for FILE in $dir; do

while [ "$(pidof $algo | wc -w)" -ge $processors ]; do
  sleep 1
done

echo ./$algo $FILE $type $filt $lr smearsumrel 1e-8 1e-8 $time $taumohc $rand $ampl $alpha $step
./$algo $FILE $type $filt $lr smearsumrel 1e-8 1e-8 $time $taumohc $rand $ampl $alpha $step &


done

#~ while [ "$(pidof $algo | wc -w)" -ge 1 ]; do
  #~ sleep 10
#~ done
