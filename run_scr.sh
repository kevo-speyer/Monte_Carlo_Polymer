#!/bin/bash
n_runs=20
i_run=1
while [ $i_run -le $n_runs ];do
    ./exe_MC_pol
    mkdir ${i_run}_run
    for fls in fort.60 fort.71 fort.73 fort.81 last_config.dat init_positions.dat ;do
        mv $fls ./${i_run}_run
    done 
    cp ${i_run}_run/last_config.dat ./init_positions.dat
    ((i_run++))
done 

for fls in fort.71 fort.60 fort.81 fort.73; do
    paste *_run/$fls  | awk '{x=0;for(i=2;i<=NF;i=i+2){x+=$i};print $1,2*x/NF}' > avg_${fls}
done

