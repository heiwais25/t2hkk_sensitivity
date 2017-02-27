#!/bin/bash

source /home/mhartz/.bashrc
cd /home/mhartz/work/t2hkk_sens/T2HKKSensitivity/Fitter

./bin/massHierarchy -c1 input_hk.card -c2 input_kd_1p5.card -o /var/tmp/mh_2det_1p5_escale_${1}.root -co osc_param.card -cs syst_param_escale.card -d $1 &> /var/tmp/mh_2det_1p5_escale_${1}.log

mv /var/tmp/mh_2det_1p5_escale_${1}.root /home/mhartz/work/t2hkk_sens/T2HKKSensitivity/Fitter
mv /var/tmp/mh_2det_1p5_escale_${1}.log /home/mhartz/work/t2hkk_sens/T2HKKSensitivity/Fitter

