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
#include "probWrapper.h"
#include "OscParams.h"
#include <omp.h>
#include <ctime>



void ParseArgs(int argc, char **argv);

char *fOutFile = NULL;
char *fOscFile = NULL;
float fBaseline;
float fDensity;

int main(int argc, char *argv[])
{

 fBaseline = 295.;
 fDensity = 2.6; 

 //Get the arguments
 ParseArgs(argc, argv);

 //Load the oscillation parameters
 OscParams *oscParams = new OscParams(fOscFile);

 //Set the oscillation probability
 ProbWrapper* prob = new ProbWrapper();
 prob->SetBaseLine(fBaseline);
 prob->SetDensity(fDensity);
 prob->SetOscillationParameters(oscParams->GetParamValue(0),oscParams->GetParamValue(1),oscParams->GetParamValue(2),oscParams->GetParamValue(3),oscParams->GetParamValue(4),oscParams->GetParamValue(5));
 

 TGraph *gr1 = new TGraph();
 TGraph *gr2 = new TGraph();
 TGraph *gr3 = new TGraph();
 TGraph *gr4 = new TGraph();
 TGraph *grmu = new TGraph();
 for(int i=0; i<2000; i++){
   float enu = 0.005+(float)i*0.005;
   grmu->SetPoint(1999-i,295.0/enu,prob->GetProbNuMuNuMu(enu));
   gr1->SetPoint(i,enu,prob->GetProbNuMuNuE(enu));
   gr2->SetPoint(i,enu,prob->GetProbNuMuBarNuEBar(enu));
   gr3->SetPoint(i,enu,prob->GetProbNuMuNuMu(enu));
   gr4->SetPoint(i,enu,prob->GetProbNuMuBarNuMuBar(enu));
 }
 TFile *fout = new TFile(fOutFile,"RECREATE");
 gr1->Write("numu_to_nue_prob");
 gr2->Write("numub_to_nueb_prob");
 gr3->Write("numu_to_numu_prob");
 gr4->Write("numub_to_numub_prob");
 grmu->Write("numu_to_numu_LoverE");
 fout->Close();
      

}   


void ParseArgs(int argc, char **argv){
  std::cout << "parse args" << std::endl;
  int nargs = 1; 
  if(argc<(nargs*2+1)){ exit(1); }
  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-co") fOscFile = argv[++i];
    if(std::string(argv[i]) == "-o") fOutFile = argv[++i];
    if(std::string(argv[i]) == "-b") fBaseline = atof(argv[++i]);
    if(std::string(argv[i]) == "-d") fDensity = atof(argv[++i]);
  } 
}

