#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=8

sigma=5
alpha=0.000625
h=1.410473959
ntraj=1000
delay=1

srvals=(0.001 0.0032 0.01 0.032 0.1 0.32 1 3.2 10)
tvals=(66 60.5 55 44 33 22 11 10 10)
len=${#srvals[@]}
echo ${srvals[*]} > "srvals.dat"
echo ${tvals[*]} > "tvals.dat"

for ((i=0; i<$len; i++)); do
 sr=${srvals[$i]}
 tmax=${tvals[$i]}
 
 echo $sr $alpha $h $sigma > "inputparameters.inp"

 echo $ntraj 4 $tmax $delay > "timestepdata.inp"
 echo 0.05 >> "timestepdata.inp"
 echo 0.02 >> "timestepdata.inp"
 echo 0.01 >> "timestepdata.inp"
 echo 0.005 >> "timestepdata.inp"

 if [ $i -lt 4 ]; then
  echo 1 0 0 > "options.inp"
 else
  echo 0 0 0 > "options.inp"
 fi

 echo "running with sr=$sr"
 time $(pwd)"/main.out"

 mv S.dat S_sr$sr.dat -f
 mv Ql.dat Ql_sr$sr.dat -f
 mv eta.dat eta_sr$sr.dat -f
 mv psi.dat psi_sr$sr.dat -f
 mv psi2.dat psi2_sr$sr.dat -f
done

