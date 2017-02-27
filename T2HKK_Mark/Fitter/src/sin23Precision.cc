#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "math.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "TF1.h"
#include "TMath.h"
#include "MakeHistograms.h"
#include "OscParams.h"
#include "SystParams.h"
#include <omp.h>


//Input function and variables
void ParseArgs(int argc, char **argv);
char *fCardFileHK = NULL;
char *fCardFileKD = NULL;
char *fCardOsc = NULL;
char *fCardSyst = NULL;
char *fOutFile = NULL;
int fHierarchy;
int fSin23;

//Function to calculate the likelihood
double CalcLogLikelihood(TH1D **data1Re, TH1D **data1Rmu, TH1D **pred1Re, TH1D **pred1Rmu, bool cut = true);

//FCN for fitter
void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  );

//Instances of class to make histograms
MakeHistograms *histsHK=NULL;
MakeHistograms *histsKD=NULL;

//Data and predicted histograms
TH1D *data1ReHK[2];
TH1D *data1RmuHK[2];
TH1D *pred1ReHK[2];
TH1D *pred1RmuHK[2];
TH1D *data1ReKD[2];
TH1D *data1RmuKD[2];
TH1D *pred1ReKD[2];
TH1D *pred1RmuKD[2];

//Oscillation probability
ProbWrapper* prob;

//Oscillation parameters
OscParams* oscParams;

//Systematic parameters
SystParams* systParams;

int main(int argc, char *argv[])
{

 //Initialize variabels
 fHierarchy = -1;
 fSin23 = -1; 

 //Get the arguments
 ParseArgs(argc, argv);

 //Set All Histograms to NULL
 for(int i=0; i<2; i++){
   data1ReHK[i] = NULL; 
   data1RmuHK[i] = NULL; 
   data1ReKD[i] = NULL; 
   data1RmuKD[i] = NULL; 
   pred1ReHK[i] = NULL;
   pred1RmuHK[i] = NULL;
   pred1ReKD[i] = NULL;
   pred1RmuKD[i] = NULL; 
 } 

 //Oscillation probability
 prob = new ProbWrapper();
 prob->SetDeltaCP(0);
 prob->SetBaseLine(1100.);
 prob->SetDensity(3.0);

 //Load the oscillation parameters
 oscParams = new OscParams(fCardOsc);

 //Load the systematic parameters 
 systParams = new SystParams(fCardSyst);

 //Class for building histograms
 histsHK = new MakeHistograms(fCardFileHK, prob, 0);
 histsHK->ApplyFluxWeights();
 if(fCardFileKD!=NULL){
   histsKD = new MakeHistograms(fCardFileKD, prob, 1);
   histsKD->ApplyFluxWeights();
 } 

 
 //Use a minuit fitter
 TFitter *minuit = new TFitter(oscParams->GetNOsc()+systParams->GetNSysts());

 //Make the output ttree
 TFile *fout = new TFile(fOutFile,"RECREATE");
 TTree results("results","CP violation significance");
 float deltacp, sin23;
 int hierarchy;
 TGraph *gr = new TGraph();
 results.Branch("DeltaChi2Scan","TGraph",&gr,32000,0);
 results.Branch("Sin23",&sin23,"Sin23/F"); // Chnaging Theta23
 results.Branch("Hierarchy",&hierarchy,"Hierarchy/I");

 //double offsets[19] = {-43,-31,-21,-17,-13,-10,-5,-3,-2,0,2,3,7,10,13,17,21,31,43}
 double offsets[33] = {-43,-32,-28,-24,-20,-18,-16,-14,-12,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,12,14,16,18,20,24,28,32,43};
 

 //Iterate through the true values of deltacp
 for(int itsin23 = 0; itsin23<101; itsin23++){
   //If a value was specified at the command line, only continue if at the value
   if(fSin23>=0 && itsin23!=fSin23) continue;
   double tsin23 = 0.5 + (double)itsin23*(0.2/100.); // Need to apply the range of sin theta 23
   sin23 = tsin23;

   //Iterate through the hierarchies
   for(int ithier = 0; ithier<2; ithier++){
     //If a value was specified at the command line, only continue if at the value
     if(fHierarchy>=0 && ithier!=fHierarchy) continue;
     hierarchy = ithier;

     //Set the oscillation parameters and settings to calculate the prediction
     prob->SetOscillationParameters((ithier==0 ? 1.0 : -1.0)*oscParams->GetParamValue(0),tsin23,oscParams->GetParamValue(2),oscParams->GetParamValue(3),oscParams->GetParamValue(4),oscParams->GetParamValue(5));

     prob->SetBaseLine(histsHK->GetBaseline());
     prob->SetDensity(histsHK->GetDensity());

     //Reset the systematics
     systParams->ResetParams();

     //Make the first detector (HK) fake data
     histsHK->ApplyOscillations();
     histsHK->GetPredictions(data1ReHK, data1RmuHK, systParams);

     //Only do the second detector if it is specified 
     if(histsKD!=NULL)
     {
       //Change the baseline and matter density
       prob->SetBaseLine(histsKD->GetBaseline());
       prob->SetDensity(histsKD->GetDensity());

       //Make the second detector (KD) fake data
       histsKD->ApplyOscillations();
       histsKD->GetPredictions(data1ReKD, data1RmuKD, systParams);
     }   

     //Reset the TGraph
     gr->Delete();
     gr = new TGraph();    
 
     //Scan the delta_cp points near the true value
     for(int isin23 = 0; isin23<33; isin23++){
       double fit_sin23 = tsin23+offsets[isin23]/400.;
       double minimum = 9e9;
       //Run the fit for both hierarchies
       for(int ihier = 0; ihier<2; ihier++){
         //Start the fit in both octants to ensure global minimum is found
         //for(int ioct = 0; ioct<2; ioct++){
      
           //Clear the fitter 
           minuit->Clear(); 

           int piter = 0;
           //Set Systematic Parameters
           for(int i=0; i<systParams->GetNSysts(); i++){
             minuit->SetParameter(piter,Form("syst%d",i),systParams->GetParamNominal(i),systParams->GetParamError(i)/10.,0,0);
             piter++;
           }

           //Set Oscillation parameters
           for(int i=0; i<oscParams->GetNOsc(); i++){
             double value = oscParams->GetParamValue(i);
             if(i==0 && ihier==1) value = -1.0*value;
             //if(i==1 && ioct==0) value = 0.40;
             //if(i==1 && ioct==1) value = 0.60;
             if(i==1) value = fit_sin23;
	    
             //Applying mass hierarchy
             if(oscParams->GetParamType(i)==0)
               minuit->SetParameter(piter,Form("osc%d",i),value,oscParams->GetParamError(i)/100.,
                                    (i==0 && ihier==1 ? -1.0*oscParams->GetParamMaximum(i) : oscParams->GetParamMinimum(i)),
                                    (i==0 && ihier==1 ? -1.0*oscParams->GetParamMinimum(i) : oscParams->GetParamMaximum(i)));
             else
               minuit->SetParameter(piter,Form("osc%d",i),value,oscParams->GetParamError(i)/100.,oscParams->GetParamMinimum(i),oscParams->GetParamMaximum(i));
             piter++;
           }
           //minuit->SetParameter(piter,Form("osc%d",piter),fit_deltacp,oscParams->GetParamError(oscParams->GetNOsc()-1)/100.0,0,0);
           minuit->FixParameter(systParams->GetNSysts());
           minuit->FixParameter(systParams->GetNSysts() + 1);

           //Set the FCN
           minuit->SetFCN(myFcn);

          
           // set print level
           double arglist[100];
           arglist[0] = 2;
           minuit->ExecuteCommand("SET PRINT",arglist,2);

           //Set the minimizer
           arglist[0] = 5000; // number of function calls
           arglist[1] = 0.001; // tolerance
           minuit->ExecuteCommand("MIGRAD",arglist,2);
           //minuit->ExecuteCommand("SIMPLEX",arglist,2);
           //minuit->ExecuteCommand("MINI",arglist,2);

           //Get the fit results
           double amin,edm,errdeff;
           int nvpar, nparx;
           int conv = minuit->GetStats(amin,edm,errdeff,nvpar,nparx);

           std::cout << "Converged: " << conv << std::endl;
           std::cout << "Minimum: " << amin << std::endl;

           //If the fit converged and the minimum is less than other starting points, save it
           if(amin<minimum && conv==3) minimum = amin;
           //if(amin<minimum && conv==0) minimum = amin;
             
         //} //octant loop
       }//hierarchy loop
       std::cout << "Minimum at " << fit_sin23 << " is " << minimum << std::endl;
       gr->SetPoint(isin23,fit_sin23,minimum);
     } //dcp loop  

  
     //Fill the ttree
     results.Fill();
   }  
 }

 //Write the ttree to the output file
 results.Write();
 fout->Close();
}

//Function to parse the input arguments
void ParseArgs(int argc, char **argv){
  std::cout << "parse args" << std::endl;
  int nargs = 1; 
  if(argc<(nargs*2+1)){ exit(1); }
  for(int i = 1; i < argc; i++){
    if(std::string(argv[i]) == "-c1") fCardFileHK = argv[++i];
    if(std::string(argv[i]) == "-c2") fCardFileKD = argv[++i];
    if(std::string(argv[i]) == "-co") fCardOsc = argv[++i];
    if(std::string(argv[i]) == "-cs") fCardSyst = argv[++i];
    if(std::string(argv[i]) == "-o") fOutFile = argv[++i];
    if(std::string(argv[i]) == "-h") fHierarchy = atoi(argv[++i]);
    if(std::string(argv[i]) == "-s") fSin23 = atoi(argv[++i]);
  } 
}

//Function to calculate the likelihood
double CalcLogLikelihood(TH1D **data1Re, TH1D **data1Rmu, TH1D **pred1Re, TH1D **pred1Rmu, bool cut){

  double ll = 0.;
  for(int i=0; i<2; i++)
    for(int j=1; j<data1Re[i]->GetNbinsX(); j++){
      double erec = data1Re[i]->GetXaxis()->GetBinCenter(j);
      if(erec<0.1) continue;
      if(erec>1.2 && cut) continue;
      if(pred1Re[i]->GetBinContent(j)>1.0e-9 && data1Re[i]->GetBinContent(j)>1.0e-9)
        ll += 2.0*(pred1Re[i]->GetBinContent(j)-data1Re[i]->GetBinContent(j)+
              data1Re[i]->GetBinContent(j)*TMath::Log( data1Re[i]->GetBinContent(j)/pred1Re[i]->GetBinContent(j)) );
      else if(data1Re[i]->GetBinContent(j)<=1.0e-9)
        ll += 2.0*pred1Re[i]->GetBinContent(j);
    } 

  for(int i=0; i<2; i++)
    for(int j=1; j<data1Rmu[i]->GetNbinsX(); j++){
      double erec = data1Rmu[i]->GetXaxis()->GetBinCenter(j);
      if(erec<0.2) continue;
      if(pred1Rmu[i]->GetBinContent(j)>1.0e-9 && data1Rmu[i]->GetBinContent(j)>1.0e-9)
        ll += 2.0*(pred1Rmu[i]->GetBinContent(j)-data1Rmu[i]->GetBinContent(j)+
              data1Rmu[i]->GetBinContent(j)*TMath::Log( data1Rmu[i]->GetBinContent(j)/pred1Rmu[i]->GetBinContent(j)) );
      else if(data1Rmu[i]->GetBinContent(j)<=1.0e-9)
        ll += 2.0*pred1Rmu[i]->GetBinContent(j);
    } 

  return ll;
}

//FCN function for the minimizer
void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  ){


  double matterDensity[2];
  matterDensity[0] = histsHK->GetDensity();
  if(histsKD!=NULL) matterDensity[1] = histsKD->GetDensity();

  //Get the updated systematic parameter values
  for(int i=0; i<systParams->GetNSysts(); i++){
    systParams->SetParamValue(i,p[i]);
    if(systParams->GetParamType(i)==MATTERDENSITY){
      if(systParams->DoDetector(1,i)) matterDensity[0] = p[i];
      if(systParams->DoDetector(2,i)) matterDensity[1] = p[i];
    }
  }
  int nSyst = systParams->GetNSysts();

  //Update the oscillation parameter values
  prob->SetOscillationParameters(p[nSyst],p[nSyst+1],p[nSyst+2],p[nSyst+3],p[nSyst+4],p[nSyst+5]);
  //Calculate the first detector prediction
  prob->SetBaseLine(histsHK->GetBaseline());
  prob->SetDensity(matterDensity[0]);
  histsHK->ApplyOscillations();
  histsHK->GetPredictions(pred1ReHK, pred1RmuHK,systParams);
  if(histsKD!=NULL){
    //Calculate the secon detector prediction
    prob->SetBaseLine(histsKD->GetBaseline());
    prob->SetDensity(matterDensity[1]);
    histsKD->ApplyOscillations();
    histsKD->GetPredictions(pred1ReKD, pred1RmuKD,systParams);
  }

  //Get the data chi2
  double chisq = 0.;
  chisq += CalcLogLikelihood(data1ReHK,data1RmuHK,pred1ReHK, pred1RmuHK, (histsHK->GetBaseline()<300.) );
  if(histsKD!=NULL) chisq += CalcLogLikelihood(data1ReKD,data1RmuKD,pred1ReKD, pred1RmuKD, (histsKD->GetBaseline()<300.) );

  //Get the oscillation parameter penalty term
  for(int i=0; i<oscParams->GetNOsc(); i++)
    if(oscParams->GetParamType(i)==1)
      chisq += pow( (p[i+nSyst]-oscParams->GetParamValue(i))/oscParams->GetParamError(i), 2);

  //Get the systematic parameter penalty term
  chisq += systParams->GetChi2();

  fval = chisq;

}


