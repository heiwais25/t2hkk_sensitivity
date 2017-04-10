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
int fDeltacp;

int czNbin = 20;
int azNbin = 12;


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

// Data and predicted histogram for carrying total value around all azimuth angle
TH1D *data1ReTotalHK[2];
TH1D *data1RmuTotalHK[2];
TH1D *pred1ReTotalHK[2];
TH1D *pred1RmuTotalHK[2];
TH1D *data1ReTotalKD[2];
TH1D *data1RmuTotalKD[2];
TH1D *pred1ReTotalKD[2];
TH1D *pred1RmuTotalKD[2];

//Data and predicted histograms
TH2D *atmData1ReHK;
TH2D *atmData1RmuHK;
TH2D *atmPred1ReHK;
TH2D *atmPred1RmuHK;
TH2D *atmData1ReKD;
TH2D *atmData1RmuKD;
TH2D *atmPred1ReKD;
TH2D *atmPred1RmuKD;

TH2D *raw1ReEvent[2][6][3];
TH2D *raw1RmuEvent[2][6][3];


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
    fDeltacp = -1; 

    //Get the arguments
    ParseArgs(argc, argv);

    //Set All Histograms to NULL
    for(int i=0; i<2; i++)
    {
        data1ReHK[i] = NULL; 
        data1RmuHK[i] = NULL; 
        data1ReKD[i] = NULL; 
        data1RmuKD[i] = NULL; 
        pred1ReHK[i] = NULL;
        pred1RmuHK[i] = NULL;
        pred1ReKD[i] = NULL;
        pred1RmuKD[i] = NULL; 

        data1ReTotalHK[i] = NULL; 
        data1RmuTotalHK[i] = NULL; 
        data1ReTotalKD[i] = NULL; 
        data1RmuTotalKD[i] = NULL; 
        pred1ReTotalHK[i] = NULL;
        pred1RmuTotalHK[i] = NULL;
        pred1ReTotalKD[i] = NULL;
        pred1RmuTotalKD[i] = NULL; 

        //Data and predicted histograms

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

 
    //Use a minuit fitter
    TFitter *minuit = new TFitter(oscParams->GetNOsc()+systParams->GetNSysts());

    //Make the output ttree
    TFile *fout = new TFile(fOutFile,"RECREATE");
    TTree results("results","CP violation significance");
    float min_known, min_unknown, deltacp;
    int hierarchy;
    results.Branch("MinKnown",&min_known,"MinKnown/F");
    results.Branch("MinUnknown",&min_unknown,"MinUnknown/F");
    results.Branch("DeltaCP",&deltacp,"DeltaCP/F");
    results.Branch("Hierarchy",&hierarchy,"Hierarchy/I");


    double czStart  = 1.0;
    double czEnd    = -1.0;
    double czStep   = (czStart - czEnd) / (double)czNbin;
    double czBin[czNbin + 1];

    for(int i = 0; i <= czNbin; i++)
    {
        czBin[i] = czEnd + 0.05 + czStep * (double)i;
    }

    double azStep = (2.0 * M_PI)/(double)azNbin;

    //Class for building histograms
    histsHK = new MakeHistograms(fCardFileHK, prob, 0);
    histsHK -> SetAngularBin(czStep, azStep);
    // histsHK->ApplyFluxWeights();
    if(fCardFileKD!=NULL)
    {
        histsKD = new MakeHistograms(fCardFileKD, prob, 1);
        histsKD -> SetAngularBin(czStep, azStep);
        // histsKD->ApplyFluxWeights();
    }  
    
    raw1ReEvent[0][0][0] = histsHK -> Getraw1Re(0, 0, 0);
    raw1RmuEvent[0][0][0] = histsHK -> Getraw1Rmu(0, 0, 0);

    // Setting basic bin structure of histogram
    atmData1ReHK = new TH2D("atmData1ReHK", "", czNbin - 1, czBin, raw1ReEvent[0][0][0]->GetYaxis()->GetNbins() - 1, raw1ReEvent[0][0][0]->GetYaxis()->GetXbins()->GetArray());
    atmData1RmuHK = new TH2D("atmData1RmuHK", "", czNbin - 1, czBin, raw1RmuEvent[0][0][0]->GetYaxis()->GetNbins() - 1, raw1RmuEvent[0][0][0]->GetYaxis()->GetXbins()->GetArray());
    atmPred1ReHK = new TH2D("atmPred1ReHK", "", czNbin - 1, czBin, raw1ReEvent[0][0][0]->GetYaxis()->GetNbins() - 1, raw1ReEvent[0][0][0]->GetYaxis()->GetXbins()->GetArray());
    atmPred1RmuHK = new TH2D("atmPred1RmuHK", "", czNbin - 1, czBin, raw1RmuEvent[0][0][0]->GetYaxis()->GetNbins() - 1, raw1RmuEvent[0][0][0]->GetYaxis()->GetXbins()->GetArray());

    atmData1ReKD = new TH2D("atmData1ReKD", "", czNbin - 1, czBin, raw1ReEvent[0][0][0]->GetYaxis()->GetNbins() - 1, raw1ReEvent[0][0][0]->GetYaxis()->GetXbins()->GetArray());
    atmData1RmuKD = new TH2D("atmData1RmuKD", "", czNbin - 1, czBin, raw1RmuEvent[0][0][0]->GetYaxis()->GetNbins() - 1, raw1RmuEvent[0][0][0]->GetYaxis()->GetXbins()->GetArray());
    atmPred1ReKD = new TH2D("atmPred1ReKD", "", czNbin - 1, czBin, raw1ReEvent[0][0][0]->GetYaxis()->GetNbins() - 1, raw1ReEvent[0][0][0]->GetYaxis()->GetXbins()->GetArray());
    atmPred1RmuKD = new TH2D("atmPred1RmuKD", "", czNbin - 1, czBin, raw1RmuEvent[0][0][0]->GetYaxis()->GetNbins() - 1, raw1RmuEvent[0][0][0]->GetYaxis()->GetXbins()->GetArray());


    //Iterate through the true values of deltacp
    for(int itdcp = 0; itdcp<101; itdcp++)
    {
        //If a value was specified at the command line, only continue if at the value
        if(fDeltacp>=0 && itdcp!=fDeltacp) continue;
        double tdcp = (double)itdcp*(3.14159*2.0/100.);
        deltacp = tdcp;

        //Iterate through the hierarchies 
        for(int ithier = 0; ithier<2; ithier++)
        {
            //If a value was specified at the command line, only continue if at the value
            if(fHierarchy>=0 && ithier!=fHierarchy) continue;
            hierarchy = ithier;

            //Set the oscillation parameters and settings to calculate the prediction
            prob->SetOscillationParameters((ithier==0 ? 1.0 : -1.0)*oscParams->GetParamValue(0),oscParams->GetParamValue(1),oscParams->GetParamValue(2),oscParams->GetParamValue(3),oscParams->GetParamValue(4),tdcp);
            // prob->SetBaseLine(histsHK->GetBaseline());
            // prob->SetDensity(histsHK->GetDensity());

            //Reset the systematics
            systParams->ResetParams();
            //Fill the 2D histogram in cosine bin and energy bin
            for(int cz = 0; cz <= czNbin; cz++)
            {

                //Summation over the azimuth angle 
                for(int az = 0; az < azNbin; az++)
                {
                    // applying anglue flux 
                    histsHK->ApplyFluxWeights(czBin[cz], az);
                    if(fCardFileKD!=NULL)
                    {
                        histsKD->ApplyFluxWeights(czBin[cz], az);
                    }  


                    //Make the first detector (HK) fake data
                    histsHK->ApplyOscillations(czBin[cz]);
                    histsHK->GetPredictions(data1ReHK, data1RmuHK, systParams);

                    //Only do the second detector if it is specified 
                    if(histsKD!=NULL)
                    {
                        //Change the baseline and matter density
                        // std::cout << "Density: " << histsKD->GetDensity() << std::endl;
                        // prob->SetBaseLine(histsKD->GetBaseline());
                        // prob->SetDensity(histsKD->GetDensity());

                        //Make the second detector (KD) fake data
                        histsKD->ApplyOscillations(czBin[cz]);
                        histsKD->GetPredictions(data1ReKD, data1RmuKD, systParams);
                    }   

                    // first azimuth angle of each cosine zeith
                    if(az == 0)
                    {
                        if(data1ReTotalHK[0]!=NULL) data1ReTotalHK[0]->Delete();
                        if(data1ReTotalHK[1]!=NULL) data1ReTotalHK[1]->Delete();
                        if(data1RmuTotalHK[0]!=NULL) data1RmuTotalHK[0]->Delete();
                        if(data1RmuTotalHK[1]!=NULL) data1RmuTotalHK[1]->Delete();
                        data1ReTotalHK[0] = (TH1D *)data1ReHK[0]->Clone();
                        data1ReTotalHK[1] = (TH1D *)data1ReHK[1]->Clone();
                        data1RmuTotalHK[0] = (TH1D *)data1RmuHK[0]->Clone();
                        data1RmuTotalHK[1] = (TH1D *)data1RmuHK[1]->Clone();
                        if(histsKD!=NULL)
                        {
                            if(data1ReTotalKD[0]!=NULL) data1ReTotalKD[0]->Delete();
                            if(data1ReTotalKD[1]!=NULL) data1ReTotalKD[1]->Delete();
                            if(data1RmuTotalKD[0]!=NULL) data1RmuTotalKD[0]->Delete();
                            if(data1RmuTotalKD[1]!=NULL) data1RmuTotalKD[1]->Delete();
                            data1ReTotalKD[0] = (TH1D *)data1ReKD[0]->Clone();
                            data1ReTotalKD[1] = (TH1D *)data1ReKD[1]->Clone();
                            data1RmuTotalKD[0] = (TH1D *)data1RmuKD[0]->Clone();
                            data1RmuTotalKD[1] = (TH1D *)data1RmuKD[1]->Clone();
                        }
                    }
                    else
                    {
                        data1ReTotalHK[0]->Add(data1ReHK[0], 1.0);
                        data1ReTotalHK[1]->Add(data1ReHK[1], 1.0);
                        data1RmuTotalHK[0]->Add(data1RmuKD[0], 1.0);
                        data1RmuTotalHK[1]->Add(data1RmuKD[1], 1.0);
                        if(histsKD!=NULL)
                        {
                            data1ReTotalKD[0]->Add(data1ReHK[0], 1.0);
                            data1ReTotalKD[1]->Add(data1ReHK[1], 1.0);
                            data1RmuTotalKD[0]->Add(data1RmuKD[0], 1.0);
                            data1RmuTotalKD[1]->Add(data1RmuKD[1], 1.0);
                        }
                    }
                }
                for (int i = 1; i <= data1ReTotalHK[0]->GetXaxis()->GetNbins(); i++)
                {
                    atmData1ReHK    -> Fill(czBin[cz], data1ReTotalHK[0]->GetXaxis()->GetBinCenter(i), data1ReTotalHK[0]->GetBinContent(i));
                    atmData1RmuHK   -> Fill(czBin[cz], data1RmuTotalKD[0]->GetXaxis()->GetBinCenter(i), data1RmuTotalKD[0]->GetBinContent(i));
                    if(histsKD!=NULL)
                    {
                        atmData1ReKD    -> Fill(czBin[cz], data1ReTotalHK[0]->GetXaxis()->GetBinCenter(i), data1ReTotalHK[0]->GetBinContent(i));
                        atmData1RmuKD   -> Fill(czBin[cz], data1RmuTotalKD[0]->GetXaxis()->GetBinCenter(i), data1RmuTotalKD[0]->GetBinContent(i));
                    }
                }
            } // End of filling data by bins of cosine zenith angle and energy bins

            //To save the minima for both hierachies
            double min_fixed = 9e9;
            double min_free = 9e9;

            //Run the fit for both hierarchies
            for(int ihier = 0; ihier<2; ihier++)
            {
                //Start the fit in both octants to ensure global minimum is found
                for(int ioct = 0; ioct<2; ioct++)
                {
                    //Run the fit for both cp conserving states
                    for(int idcp = 0; idcp<2; idcp++)
                    {

                        //Clear the fitter 
                        minuit->Clear(); 

                        int piter = 0;
                        //Set Systematic Parameters
                        for(int i=0; i<systParams->GetNSysts(); i++)
                        {
                            minuit->SetParameter(piter,Form("syst%d",i),systParams->GetParamNominal(i),systParams->GetParamError(i)/10.,0,0);
                            piter++;
                        }

                        //Set Oscillation parameters
                        for(int i=0; i<oscParams->GetNOsc()-1; i++)
                        {
                            double value = oscParams->GetParamValue(i);
                            if(i==0 && ihier==1) value = -1.0*value;
                            if(i==1 && ioct==0) value = 0.40;
                            if(i==1 && ioct==1) value = 0.60;
                            if(oscParams->GetParamType(i)==0)
                            minuit->SetParameter(piter,Form("osc%d",i),value,oscParams->GetParamError(i)/100.,
                            (i==0 && ihier==1 ? -1.0*oscParams->GetParamMaximum(i) : oscParams->GetParamMinimum(i)),
                            (i==0 && ihier==1 ? -1.0*oscParams->GetParamMinimum(i) : oscParams->GetParamMaximum(i)));
                            else
                            minuit->SetParameter(piter,Form("osc%d",i),value,oscParams->GetParamError(i)/100.,oscParams->GetParamMinimum(i),oscParams->GetParamMaximum(i));
                            piter++;
                        }
                        minuit->SetParameter(piter,Form("osc%d",oscParams->GetNOsc()-1),3.14159*(float)idcp,oscParams->GetParamError(oscParams->GetNOsc()-1)/100.0,
                                        oscParams->GetParamMinimum(oscParams->GetNOsc()-1),oscParams->GetParamMaximum(oscParams->GetNOsc()-1));
                        minuit->FixParameter(piter);

                        //Set the FCN
                        minuit->SetFCN(myFcn);


                        // set print level
                        double arglist[100];
                        arglist[0] = 2;
                        minuit->ExecuteCommand("SET PRINT",arglist,2);

                        //Set the minimizer
                        arglist[0] = 50000; // number of function calls
                        arglist[1] = 0.001; // tolerance
                        minuit->ExecuteCommand("MIGRAD",arglist,2);
                        //minuit->ExecuteCommand("MINI",arglist,2);

                        //Get the fit results
                        double amin,edm,errdeff;
                        int nvpar, nparx;
                        int conv = minuit->GetStats(amin,edm,errdeff,nvpar,nparx);

                        //If the fit converged and the minimum is less than other starting points, save it
                        if(amin<min_free && conv>0) min_free = amin;
                        if(amin<min_fixed && (ihier==ithier) && conv>0) min_fixed = amin;
                         
                    } //dcp loop  
                } //octant loop
            }//hierarchy loop

            //Save the minima and deltachi2
            min_known = min_fixed;
            min_unknown = min_free;

            //Fill the ttree
            results.Fill();
        }  // ithierarchy loop
    } // itdcp loop
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
    if(std::string(argv[i]) == "-d") fDeltacp = atoi(argv[++i]);
  } 
}

//Function to calculate the likelihood
double CalcLogLikelihood(TH1D **data1Re, TH1D **data1Rmu, TH1D **pred1Re, TH1D **pred1Rmu, bool cut){

  double ll = 0.;
  for(int i=0; i<2; i++)
    for(int j=1; j<data1Re[i]->GetNbinsX(); j++){
      double erec = data1Re[i]->GetXaxis()->GetBinCenter(j);
      if(erec>1.2 && cut) continue;
      if(erec<0.1) continue;
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

double CalcLogLikelihood2D(TH2D *data1ReTotal, TH2D *data1RmuTotal, TH2D *pred1ReTotal, TH2D *pred1RmuTotal, bool cut)
{

    double ll = 0.;
    for(int j=1; j<data1ReTotal->GetXaxis()->GetNbins(); j++) // Cosine axis
    {
        for(int k=1; k<data1ReTotal->GetYaxis()->GetNbins(); k++) // Energy axis
        {
            double erec = data1ReTotal->GetYaxis()->GetBinCenter(j);
            if(erec>1.2 && cut) continue;
            if(erec<0.1) continue;
            if(pred1ReTotal->GetBinContent(j,k)>1.0e-9 && data1ReTotal->GetBinContent(j,k)>1.0e-9)
                ll += 2.0*(pred1ReTotal->GetBinContent(j,k)-data1ReTotal->GetBinContent(j,k)+
                    data1ReTotal->GetBinContent(j,k)*TMath::Log( data1ReTotal->GetBinContent(j,k)/pred1ReTotal->GetBinContent(j,k)) );
            else if(data1ReTotal->GetBinContent(j,k)<=1.0e-9)
                ll += 2.0*pred1ReTotal->GetBinContent(j,k);
        }
    } 

    for(int j=1; j<data1RmuTotal->GetXaxis()->GetNbins(); j++) // Cosine axis
    {
        for(int k=1; k<data1RmuTotal->GetYaxis()->GetNbins(); k++) // Energy axis
        {
            double erec = data1RmuTotal->GetYaxis()->GetBinCenter(j);
            if(erec<0.2) continue;
            if(pred1RmuTotal->GetBinContent(j,k)>1.0e-9 && data1RmuTotal->GetBinContent(j,k)>1.0e-9)
                ll += 2.0*(pred1RmuTotal->GetBinContent(j,k)-data1RmuTotal->GetBinContent(j,k)+
                    data1RmuTotal->GetBinContent(j,k)*TMath::Log( data1RmuTotal->GetBinContent(j,k)/pred1RmuTotal->GetBinContent(j,k)) );
            else if(data1RmuTotal->GetBinContent(j,k)<=1.0e-9)
                ll += 2.0*pred1RmuTotal->GetBinContent(j,k);
        }
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
    // prob->SetBaseLine(histsHK->GetBaseline());
    // prob->SetDensity(matterDensity[0]);

    // Atmoshperic neutrino flux data
    double czStart  = 1.0;
    double czEnd    = -1.0;
    double czStep   = (czStart - czEnd) / (double)czNbin;
    double czBin[czNbin + 1];

    for(int i = 0; i <= czNbin; i++)
    {
        czBin[i] = czEnd + 0.05 + czStep * (double)i;
    }
    double azStep = (2.0 * M_PI)/(double)azNbin;

    for(int cz = 0; cz < czNbin; cz++)
    {
        for(int az = 0; az < azNbin; az++)
        {
            histsHK->ApplyFluxWeights(czBin[cz], az);
            histsHK->ApplyOscillations(czBin[cz]);
            histsHK->GetPredictions(pred1ReHK, pred1RmuHK,systParams);
            if(histsKD!=NULL)
            {
                //Calculate the secon detector prediction
                // prob->SetBaseLine(histsKD->GetBaseline());
                // prob->SetDensity(matterDensity[1]);
                histsKD->ApplyFluxWeights(czBin[cz], az);
                histsKD->ApplyOscillations(czBin[cz]);
                histsKD->ApplyOscillations();
                histsKD->GetPredictions(pred1ReKD, pred1RmuKD,systParams);
            }        

            // first azimuth angle of each cosine zeith
            if(az == 0)
            {
                if(pred1ReTotalHK[0]!=NULL) pred1ReTotalHK[0]->Delete();
                if(pred1ReTotalHK[1]!=NULL) pred1ReTotalHK[1]->Delete();
                if(pred1RmuTotalHK[0]!=NULL) pred1RmuTotalHK[0]->Delete();
                if(pred1RmuTotalHK[1]!=NULL) pred1RmuTotalHK[1]->Delete();
                pred1ReTotalHK[0] = (TH1D *)pred1ReHK[0]->Clone();
                pred1ReTotalHK[1] = (TH1D *)pred1ReHK[1]->Clone();
                pred1RmuTotalHK[0] = (TH1D *)pred1RmuHK[0]->Clone();
                pred1RmuTotalHK[1] = (TH1D *)pred1RmuHK[1]->Clone();
                if(histsKD!=NULL)
                {
                    if(pred1ReTotalKD[0]!=NULL) pred1ReTotalKD[0]->Delete();
                    if(pred1ReTotalKD[1]!=NULL) pred1ReTotalKD[1]->Delete();
                    if(pred1RmuTotalKD[0]!=NULL) pred1RmuTotalKD[0]->Delete();
                    if(pred1RmuTotalKD[1]!=NULL) pred1RmuTotalKD[1]->Delete();
                    pred1ReTotalKD[0] = (TH1D *)pred1ReKD[0]->Clone();
                    pred1ReTotalKD[1] = (TH1D *)pred1ReKD[1]->Clone();
                    pred1RmuTotalKD[0] = (TH1D *)pred1RmuKD[0]->Clone();
                    pred1RmuTotalKD[1] = (TH1D *)pred1RmuKD[1]->Clone();
                }
            }
            else
            {
                pred1ReTotalHK[0]->Add(pred1ReHK[0], 1.0);
                pred1ReTotalHK[1]->Add(pred1ReHK[1], 1.0);
                pred1RmuTotalHK[0]->Add(pred1RmuHK[0], 1.0);
                pred1RmuTotalHK[1]->Add(pred1RmuHK[1], 1.0);
                if(histsKD!=NULL)
                {
                    pred1ReTotalKD[0]->Add(pred1ReKD[0], 1.0);
                    pred1ReTotalKD[1]->Add(pred1ReKD[1], 1.0);
                    pred1RmuTotalKD[0]->Add(pred1RmuKD[0], 1.0);
                    pred1RmuTotalKD[1]->Add(pred1RmuKD[1], 1.0);
                }
            }

        } // END of azimuth angle
        
        for (int i = 1; i <= pred1ReTotalHK[0]->GetXaxis()->GetNbins(); i++)
        {
            atmPred1ReHK    -> Fill(czBin[cz], pred1ReTotalHK[0]->GetXaxis()->GetBinCenter(i), pred1ReTotalHK[0]->GetBinContent(i));
            atmPred1RmuHK   -> Fill(czBin[cz], pred1RmuTotalKD[0]->GetXaxis()->GetBinCenter(i), pred1RmuTotalKD[0]->GetBinContent(i));
            if(histsKD!=NULL)
            {
                atmPred1ReKD    -> Fill(czBin[cz], pred1ReTotalHK[0]->GetXaxis()->GetBinCenter(i), pred1ReTotalHK[0]->GetBinContent(i));
                atmPred1RmuKD   -> Fill(czBin[cz], pred1RmuTotalKD[0]->GetXaxis()->GetBinCenter(i), pred1RmuTotalKD[0]->GetBinContent(i));
            }
        }
    }   // END of cosine zenith angle

    //Get the data chi2
    double chisq = 0.;
    chisq += CalcLogLikelihood2D(atmData1ReHK,atmData1RmuHK,atmPred1ReHK, atmPred1RmuHK, (histsHK->GetBaseline()<300.) );
    if(histsKD!=NULL) chisq += CalcLogLikelihood2D(atmData1ReKD,atmData1RmuKD,atmPred1ReKD, atmPred1RmuKD, (histsKD->GetBaseline()<300.) );
  

    //Get the oscillation parameter penalty term
    for(int i=0; i<oscParams->GetNOsc(); i++)
        if(oscParams->GetParamType(i)==1)
            chisq += pow( (p[i+nSyst]-oscParams->GetParamValue(i))/oscParams->GetParamError(i), 2);

    //Get the systematic parameter penalty term
    chisq += systParams->GetChi2();

    fval = chisq;

}
