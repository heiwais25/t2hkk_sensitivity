#!/bin/bash

source /home/mhartz/.bashrc
cd /home/mhartz/work/t2hkktmp/Fitter

var1=$(expr $1 / 30)
var2=$(expr $1 % 30)

var2=`echo "0.35+${var2}*0.01" | bc`
var1=`echo "0.00235+${var1}*0.00001" | bc`

#-cs syst_param.card
./bin/min2D -c1 input_hk.card -co osc_param.card -p1 ${var1} -n1 dm32 -p2 ${var2} -n2 ssqth23 -fh 1 -o /var/tmp/dmsq32_ssqth32_hk_2p5_nosyst_${1}.root &> /var/tmp/dmsq32_ssqth32_hk_2p5_nosyst_${1}.log

mv /var/tmp/dmsq32_ssqth32_hk_2p5_nosyst_${1}.root /home/mhartz/work/T2HKK/Sensitivity/New/
mv /var/tmp/dmsq32_ssqth32_hk_2p5_nosyst_${1}.log /home/mhartz/work/T2HKK/Sensitivity/New/

