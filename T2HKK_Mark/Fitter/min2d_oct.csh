#!/bin/tcsh -f 
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH":/reno/home/heiwais25/T2HKKSensitivity/Fitter/lib"

set workdir = /reno/home/heiwais25/T2HKKSensitivity/Fitter 
cd $workdir
#@ one = 80/40
# 40 bin #
#@ first = $1 / 40
#@ second = $1 % 40
# 100 bin #
@ first = $1 / 100
@ second = $1 % 100
#echo "$first + $second"
#@ first = $1
#calc "2*pi"
set Dect1 = hk
#echo "$one $other"
set Dect2 = mt_bisul
# // hk // kd_1p5 // kd_2p0 // kd_2p5 // mt_bisul

# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/output 
# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/output_both_dect
# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/for_kps
#set outdir = /reno/home/sykim/study_T2HKK/T2HKK/comp_nova_min
set outdir = /reno/home/heiwais25/T2HKKSensitivity/combine_source/
# set outdir = /reno/home/sykim/study_T2HKK/T2HKK/comp_nova_max
set Type_mode = min2D_oct
 #  // cpPrecision  // cpViolation // massHierarchy // sin23Precision // sin23Violation // min2D
 #  // min2D_new    // min2D_oct
#set output = sin23_prec_mt_bisul.root

#./bin/$Type_mode -c1 input_$Dect2.card  -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param.card -cs mod_syst_param_total.card -d $1

#./bin/$Type_mode -c1 input_$Dect2.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param.card -cs syst_param_total2.card -d $1

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param_nova.card -cs syst_param_total.card -d $1

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param.card -cs syst_param_total.card -d $1

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/${Type_mode}_output_run_$1.root -co osc_param.card -cs syst_param_total.card -s $1

# change 2 oscillation parameter
#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect1/${Type_mode}_sin23_0.50_${first}_dm23_0.0024_${second}_run_gr_10000_2hk_reactor.root -co osc_param.card -cs syst_param_total.card -p1 0.0024 -p2 0.50 -n1 dm32 -n2 ssqth23 -s1 ${first} -s2 ${second}

# change 1 oscillation parameter
#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/${Type_mode}_1oct_${first}_2oct_${second}_run_10000_2hk_reactor_total.root -co osc_param.card -cs syst_param_total.card -s1 ${first} -s2 ${second}
./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/source_${first}_${second}_run_1845.root -co osc_param.card -cs syst_param_total.card -s1 ${first} -s2 ${second}


# Only one detector case
#./bin/$Type_mode -c1 input_$Dect1.card -o $outdir/$Type_mode/$Dect1/${Type_mode}_1oct_${first}_2oct_${second}_run_1839.root -co osc_param.card -cs syst_param_total_KDonly.card -s1 ${first} -s2 ${second}
#./bin/$Type_mode -c1 input_$Dect1.card -o $outdir/source_${first}_${second}_run_2009.root -co osc_param.card -cs syst_param_total_KDonly.card -s1 ${first} -s2 ${second}

#./bin/$Type_mode -c1 input_$Dect1.card -c2 input_$Dect2.card -o $outdir/$Type_mode/$Dect2/output_run_$1.root -co osc_param_nova_m.card -cs syst_param_total.card -d $1

