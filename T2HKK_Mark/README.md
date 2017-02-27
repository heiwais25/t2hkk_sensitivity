# T2HKKSensitivity
Sensitivity code for T2HKK

To compile it:
-------------

cd Prob3++ 

make

cd ../Fitter

make links 

make

The code is now compiled
Configuration card files:
----------------------------
In the fitter directory there are .card files used to configure the fits for sensitivity studies.  The input_*.card files configure a detector:

input_hk.card = Hyper-K at Tochibora configuration
input_kd_1p5.card = Korean detector at 1.5 degrees off-axis
input_kd_2p0.card = Korean detector at 2.0 degrees off-axis
input_kd_2p5.card = Korean detector at 2.0 degrees off-axis

The other two card files are osc_param.card and syst_param.card.  They configure the oscillation parameters and systematic parameters respectively.  The format for the card files is described here:

input_*.card - The first line specifies the flux weights files that is used to reweight the 295 km 2.5 degree off-axis flux.  The next 24 lines contain the POT weights and the root files that contain the event rate predictions in bins of true neutrino energy and reconstructed neutrino energy.  The histograms in each file are weighted by the numbers above to have the proper normalization. Next there are 5 numbers giving: the baseline in km, the average matter density in g/cm^3, the fiducial mass in kilotons, the neutrino mode POT in units of 1e21 and the antineutrino mode POT in units of 1e21.  To work properly, the order of all lines in this file must be maintained.

osc_param.card - Each line (with a comment # character) contains the: value, uncertainty, type, minimum value and maximum value for each parameter.  The value is the central value for the parameter, and the error should be set to any external constraint, i.e. PDG values.  The type is 0 if the parameter has a flat prior and 1 if it has a Gaussian prior (external constraint).  For parameters with flat priors, the range must also be specified with the Minimum and Maximum fields.  The oscillation parameters must be listed in the order shown in the file.

syst_param.car - Each line contains the following categories:
Value   Error   Emin  Emax  FluxMode  XsecMode  PolarityMode Detectors  Type  Correlations
Value = mean value of the parameter
Error = 1 sigma uncertainty on the parameter
Emin = equals the minimum kinematic value for which the parameter should be applied.  If Type=0, the variable is reconstructed energy, if Type=1 the variable is true neutrino energy.  Units are GeV.  This feature is not yet implemented.
Emax = equals the maximum kinematic value for which the parameter should be applied.  If Type=0, the variable is reconstructed energy, if Type=1 the variable is true neutrino energy.  Units are GeV.  This feature is not yet implemented.
FluxMode = a bit map that masks the flux channels to which the parameter should be applied.  From right to left, they are:
numu->numu, numu-bar->numu-bar, nue->nue, nue-bar->nue-bar, numu->nue, numu-bar->nue-bar.  The parameter will be applied if the corresponding bit is set to 1.
XsecMode = a bit map that masks the cross section channels to which the parameter should be applied.  From right to left, they are CCQE, CCnonQE, NC.
PolarityMode = a bit map that masks the horn polarity mode to which the parameter should be applied.  From right to left, they are neutrino mode and antineutrino mode.
Detectors = a bit map that masks detectors at which the parameter should be applied.  From right to left, they are the detector specified with the -c1 argument on the command line and the detector specified with the -c2 argument on the command line.
Type = the type of systematic.  Currently there are three: normalization in bins of reconstructed energy, normalization in bins of true neutrino energy, and an energy scale parameter.  The energy scale parameter is still under development.
Correlations = list the correlation of this parameter with all other parameters.  There should be a value for each parameter with the top-to-bottom list of parameters corresponding to the left-to-right list of correlations.

Running executables:
--------------------
After compiling, the executables in the bin directory are:
cpPrecision
cpViolation
massHierarchy
saveHistograms
The corresponding code can be found in the src directory.  Example command line executions are shown here:

cpPrecision:

./bin/cpPrecision -c1 input_hk.card -c2 input_kd_1p5.card -o output.root -co osc_param.card -cs syst_param.card -d 0 &> stdout.log

The -c1 and -c2 options are used to specify the detectors for which the sensitivities should be run.  If -c2 is not specified, the sensitivity is only evaluated with one detector.  The -co and -cs options are used to specify the oscillation parameter and systematic error files.  The -d option is used to specified which delta_cp value for which the precision should be calculated.  The value of delta_cp is 2*pi/100*<-d argument>.  The argument of -d can take values from 0 to 100.  By using the -d option it is possible to submit parallel jobs to a cluster to calculate the precision at each delta_cp value. There is also a -h option to specify true normal (0) or inverted (1) hierarchy.  If not specified, the precision for both will be calculated.  The output file contains a TTree with the name results, that has 3 branches.  DeltaChi2Scan is the profiled Delta-chi^2 curve near the true value, which can be used to calculate the precision.   DeltaCP is the true delta_cp value and Hierarchy is the true hierarchy (0=NH, 1=IH).  

cpViolation:

./bin/cpViolation -c1 input_hk.card -c2 input_kd_1p5.card -o output.root -co osc_param.card -cs syst_param.card -d 0 &> stdout.log

The arguments are essentially the same as the cpPrecision executable.  In this case the output file TTree does not contain the DeltaChi2Scan variable, but rather two variables: MinKnown and MinUnknown.  The MinKnown variable contains the minimum Delta-chi^2 among the two CP conserving states of the true known hierarchy.  The MinUnKnown contains the minimum Delta-chi^2 among the 4 CP conserving states from both hierarchies.  

massHierarchy:

./bin/massHierarchy -c1 input_hk.card -c2 input_kd_1p5.card -o output.root -co osc_param.card -cs syst_param.card -d 0 &> stdout.log

Once again, the command line arguments are the same.  There are three new variables in the output TTree.  MinNH and MinIH contain the minimum Delta-chi^2 for the fits with NH or IH fixed.  For the hierarchy fixed to the true hierarchy, this value should be 0.  It should be non-zero for the not true hierarchy.  The DeltaChi2 value is the Delta-chi^2 between the fits with the two different hierarchies.

saveHistograms:

./bin/saveHistograms -c input_hk.card -co osc_param.card

This saves a number of output root files with the predicted histograms calculated at various delta_cp and hierarchy values.  The rest of the parameters are specified through osc_param.card.
