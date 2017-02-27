#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include "math.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"
#include "OscParams.h"
#include "SystParams.h"
#include "MakeHistograms.h"
#include <omp.h>



void ParseArgs(int argc, char **argv);

void PrintMean(float *array, int nentries);

char *fCardFile = NULL;
char *fCardOsc = NULL;
char *fCardSyst = NULL;

int main(int argc, char *argv[])
{

 //Get the arguments
 ParseArgs(argc, argv);

 //Set the oscillation probability
 ProbWrapper* prob = new ProbWrapper();

 OscParams* oscParams = new OscParams(fCardOsc);
 prob->SetOscillationParameters(oscParams->GetParamValue(0),oscParams->GetParamValue(1),
                                oscParams->GetParamValue(2),oscParams->GetParamValue(3),oscParams->GetParamValue(4),oscParams->GetParamValue(5)*0.);

 SystParams *systParams = new SystParams(fCardSyst); 

 //Class for building histograms
 MakeHistograms *hists = new MakeHistograms(fCardFile, prob);
 hists->ApplyFluxWeights();
 prob->SetBaseLine(hists->GetBaseline());
 prob->SetDensity(hists->GetDensity()); 

 TH1D *pred1Re[2] = {NULL,NULL};
 TH1D *pred1Rmu[2] = {NULL,NULL};

 float int1ReFHC[10000];
 float int1ReRHC[10000];
 float int1RmuFHC[10000];
 float int1RmuRHC[10000];
 float int1ReRatio[10000];
 float erec1ReFHC[10000];
 float erec1ReRHC[10000];
 float erec1RmuFHC[10000];
 float erec1RmuRHC[10000];

 for(int i=0; i<1000; i++){

   //Throw systematic parameters
   systParams->ThrowParameters();

   //Get the matter density
   double matterDensity = hists->GetDensity();
   for(int si=0; si<systParams->GetNSysts(); si++)
     if(systParams->GetParamType(si)==MATTERDENSITY) matterDensity = systParams->GetParamValue(si);

   if( hists->GetBaseline()>300.) prob->SetDensity(matterDensity); 

   hists->ApplyOscillations();
   hists->GetPredictions(pred1Re, pred1Rmu,systParams);
   if(hists->GetBaseline()>300.){
    int1ReFHC[i] = pred1Re[0]->Integral();
    int1ReRHC[i] = pred1Re[1]->Integral();
   } else {
    int1ReFHC[i] = pred1Re[0]->Integral(1,pred1Re[0]->GetXaxis()->FindBin(1.2));
    int1ReRHC[i] = pred1Re[1]->Integral(1,pred1Re[1]->GetXaxis()->FindBin(1.2));
   }
   int1RmuFHC[i] = pred1Rmu[0]->Integral();
   int1RmuRHC[i] = pred1Rmu[1]->Integral();
   int1ReRatio[i] = int1ReFHC[i]/int1ReRHC[i];
   erec1ReFHC[i] = pred1Re[0]->GetMean();
   erec1ReRHC[i] = pred1Re[1]->GetMean();
   erec1RmuFHC[i] = pred1Rmu[0]->GetMean();
   erec1RmuRHC[i] = pred1Rmu[1]->GetMean();
   //std::cout << pred1Re[0]->Integral() << " " << pred1Re[1]->Integral() << " " << pred1Rmu[0]->Integral() << " " << pred1Rmu[1]->Integral() << std::endl;
 } 

 PrintMean(int1RmuFHC, 1000);
 PrintMean(int1RmuRHC, 1000);
 PrintMean(int1ReFHC, 1000);
 PrintMean(int1ReRHC, 1000);
 PrintMean(int1ReRatio, 1000);
 PrintMean(erec1RmuFHC, 1000);
 PrintMean(erec1RmuRHC, 1000);
 PrintMean(erec1ReFHC, 1000);
 PrintMean(erec1ReRHC, 1000);

}

void ParseArgs(int argc, char **argv){
  std::cout << "parse args" << std::endl;
  int nargs = 1; 
  if(argc<(nargs*2+1)){ exit(1); }
  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-c") fCardFile = argv[++i];
    if(std::string(argv[i]) == "-co") fCardOsc = argv[++i];
    if(std::string(argv[i]) == "-cs") fCardSyst = argv[++i];

  } 
}

void PrintMean(float *array, int nentries){

  float mean = 0.;
  for(int i=0; i<nentries; i++) mean += array[i];
  mean = mean/(float)nentries;
  float rms = 0.;
  for(int i=0; i<nentries; i++) rms += pow(array[i]-mean,2);
  rms = sqrt(rms/(float)nentries);

  std::cout << "Mean: " << mean << "   RMS: " << rms/mean*100. << std::endl;
  return;

} 
