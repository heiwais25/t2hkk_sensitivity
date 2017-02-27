#!/bin/bash

for n in {0..10}
  do
    cp run_sin23.sh run_sin23_$n.sh
    sed -i 's/{$1}/'$n'/g' run_sin23_$n.sh
    qsub -q reno run_sin23_$n.sh
    #sleep 1
    #rm -f run_sin23_$n.sh
  done

