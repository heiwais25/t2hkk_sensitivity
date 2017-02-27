#!/bin/bash

source /home/mhartz/.bashrc
cd /home/mhartz/work/t2hkk_sens/T2HKKSensitivity/Fitter/

#c2 input_kd_1p5.card
./bin/cpViolation -c1 input_hk.card -c2 input_hk.card -o /var/tmp/cpv_2det_hk_total_${1}.root -co osc_param.card -cs syst_param_total_hk.card -d $1 &> /var/tmp/cpv_2det_hk_total_${1}.log

mv /var/tmp/cpv_2det_hk_total_${1}.root /home/mhartz/work/t2hkk_sens/T2HKKSensitivity/Fitter/
mv /var/tmp/cpv_2det_hk_total_${1}.log /home/mhartz/work/t2hkk_sens/T2HKKSensitivity/Fitter/

