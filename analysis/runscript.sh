#!/bin/bash
#

mkdir -p ../data

for beta in 0.94 0.95 0.96 0.97  0.98 0.99 1.0 1.01 1.02;
do
    echo ${beta}
    if [[ ! -f ../data/beta${beta}.txt ]];
    then
        ../code/exec ${beta} > ../data/beta${beta}.txt
    fi
done

