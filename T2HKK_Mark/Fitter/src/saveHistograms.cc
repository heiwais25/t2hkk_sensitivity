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

 SystParams *systParams = new SystParams(fCardSyst); 

 //Class for building histograms
 MakeHistograms *hists = new MakeHistograms(fCardFile, prob);
 hists->ApplyFluxWeights();
 prob->SetBaseLine(hists->GetBaseline());
 prob->SetDensity(hists->GetDensity()); 

 char hier[2][20] = {"nh","ih"};
 //int dcp[8] = {0,45,90,135,180,225,270,315};
 int dcp[4] = {0,90,180,270};

 TH1D *pred1Re[2] = {NULL,NULL};
 TH1D *pred1Rmu[2] = {NULL,NULL};

 for(int i=0; i<2; i++)
   for(int j=0; j<4; j++){
     prob->SetOscillationParameters((i==0 ? 1.0 : -1.0)*oscParams->GetParamValue(0),oscParams->GetParamValue(1),
                                    oscParams->GetParamValue(2),oscParams->GetParamValue(3),oscParams->GetParamValue(4),(float)dcp[j]*TMath::Pi()/180.0);
     hists->BuildHistograms();
     hists->SaveToFile(Form("histograms_dcp%d_%s.root",dcp[j],hier[i]));
     hists->ApplyOscillations();
     hists->GetPredictions(pred1Re, pred1Rmu,systParams);
     TFile *fout = new TFile(Form("total_hists_dcp%d_%s.root",dcp[j],hier[i]),"RECREATE");
     pred1Re[0]->Write("histNuMode1Re");
     pred1Re[1]->Write("histANuMode1Re");
     pred1Rmu[0]->Write("histNuMode1Rmu");
     pred1Rmu[1]->Write("histANuMode1Rmu");
     fout->Close();
   } 

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

