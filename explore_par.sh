#!/bin/bash

for a in 0.01 0.05 0.1 0.5 1 2 5 20 100 500; do 
    awk -v a_par=$a '/a_box/{print a_par,"   ",$2}!/a_box/' input.dat > input.dat_tmp; 
    mv input.dat_tmp input.dat; ./exe_MC_pol ; 
    for fls in fort.60 fort.70 fort.71; do
        mv $fls ${fls}_a$a
    done
done

xmgrace $(for a in 0.01 0.05 0.1 0.5 1 2 5 20 100 500; do echo fort.60_a$a;done) &
