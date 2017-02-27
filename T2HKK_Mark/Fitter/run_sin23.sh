#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/reno/home/heiwais25/T2HKKSensitivity/Fitter/lib
source /reno/home/heiwais25/.bashrc
cd /reno/home/heiwais25/T2HKKSensitivity/Fitter/

source_file1=/reno/home/heiwais25/T2HKKSensitivity/Fitter/input_hk.card
source_file2=/reno/home/heiwais25/T2HKKSensitivity/Fitter/input_mt_bisul.card
output=/var/tmp/sin23_prec_2det_bisul_{$1}.root
log=/var/tmp/sin23_prec_2det_bisul_{$1}.log
osc_card=/reno/home/heiwais25/T2HKKSensitivity/Fitter/osc_param.card
sys_card=/reno/home/heiwais25/T2HKKSensitivity/Fitter/syst_param_total.card
output_dir=/reno/home/heiwais25/T2HKKSensitivity/combine_source/sin23pp/2det_bisul

#c2 input_kd_1p5.card
./bin/sin23Precision -c1 ${source_file1} -c2 ${source_file2} -o ${output} -co ${osc_card} -cs ${sys_card} -s {$1} 
 
mv ${output} ${output_dir}
#mv ${log} ${output_dir}

