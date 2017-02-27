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
//dm32, ssqth23, dm21, ssqth13, ssqth12, dcp
int fFixHierarchy;
int fSet1;
int fSet2;
double DensityKD;

// Selected value
float fParam1, fParam2;


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
    fFixHierarchy = 1;
    fSet1 = -1;
    fSet2 = -1;
    DensityKD = 0;

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
    } 

    //Oscillation probability
    prob = new ProbWrapper();
    prob->SetDeltaCP(0);
    prob->SetBaseLine(1100.);
    prob->SetDensity(3.0);

    oscParams = new OscParams(fCardOsc);

    systParams = new SystParams(fCardSyst);
    //Class for building histograms
    histsHK = new MakeHistograms(fCardFileHK, prob, 0);
    histsHK->ApplyFluxWeights();
    if(fCardFileKD!=NULL)
    {
        histsKD = new MakeHistograms(fCardFileKD, prob, 1);
        histsKD->ApplyFluxWeights();
        DensityKD = histsKD->GetDensity();

        std::cout << "Two detector is included" << std::endl;
    }

    // Starting points for profiled oscillation parameters
    int npoints[6] = {1,1,1,1,1,1};
    float startPoints[6][16]; 
    // Set initial osciallation parameter to each startPoints
    startPoints[0][0] = oscParams->GetParamValue(0);
    startPoints[1][0] = oscParams->GetParamValue(1);
    startPoints[2][0] = oscParams->GetParamValue(2);
    startPoints[3][0] = oscParams->GetParamValue(3);
    startPoints[4][0] = oscParams->GetParamValue(4);
    startPoints[5][0] = oscParams->GetParamValue(5);

    //save the information of the fixed parameters
    bool fixedParam[6] = {false, false, false, false, false, false};
    float StandardValue = 0;
    float CompareValue  = 0;

    // true value in first octant
    StandardValue = 0.35 + (double)fSet1 / 100. * 0.3; // 100bin (-+ 0.2)
    fixedParam[1] = true;
    npoints[1] = 1;
    std::cout << "first octant fixed value  = " << StandardValue << std::endl; 

    // Comparing data
    if(StandardValue < 0.5)
    {
        CompareValue = 1.0 - (double)fSet2 / 100. * 0.5;// 100bin (-+ 0.2)
    }
    else if(StandardValue > 0.5)
    {
        CompareValue = (double)fSet2 / 100. * 0.5;
    }
    else if(StandardValue == 0.5)
    {
        CompareValue = 0.5;
    }
    startPoints[1][0] = CompareValue;
    std::cout << "second octant data value  = " << CompareValue << std::endl; 

    /* Setting hierarchy */
    if(fFixHierarchy == 0) 
    {
        	startPoints[0][0] = -1.0*startPoints[0][0];
    }

    std::cout << 
        "Oscialltion Parameter in calculating data " << std::endl <<
        "dm32  : " << startPoints[0][0] << std::endl <<
        "sin23 : " << startPoints[1][0] << std::endl <<
        "dm21  : " << oscParams->GetParamValue(2) << std::endl <<
        "sin13 : " << oscParams->GetParamValue(3) << std::endl <<
        "sin12 : " << oscParams->GetParamValue(4) << std::endl <<
        "dcp   : " << oscParams->GetParamValue(5) << std::endl;

    /*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Calculating fake experimental data(data1ReHK, data1RmuHK, data1ReKD, data1RmuKD)
        The parameters which applied in this calculation is set in the oscParams file(default)
    ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
    //Set the oscillation parameters and settings to calculate the prediction
    prob->SetOscillationParameters(oscParams->GetParamValue(0),startPoints[1][0],oscParams->GetParamValue(2),
     oscParams->GetParamValue(3),oscParams->GetParamValue(4),oscParams->GetParamValue(5));
    prob->SetBaseLine(histsHK->GetBaseline());
    prob->SetDensity(histsHK->GetDensity());
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
    /*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        End of Calculating fake experimental data
    ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/


    //Use a minuit fitter
    TFitter *minuit = new TFitter(oscParams->GetNOsc()+systParams->GetNSysts());

    //Make the output ttree
    TFile *fout = new TFile(fOutFile,"RECREATE");
    TTree results("results","Minimum with 2 parameters fixed");
    float minimum =  9e9;
    float p1value = 0, p2value = 0; // p1value(standard), p2value(compare)
    int CurrentRunNumber;

    TGraph *gr = new TGraph();
    TGraph *tr = new TGraph();

    results.Branch("Minimum",&minimum,"Minimum/F");
    results.Branch("MinimumP1Value","TGraph",&gr,32000,0);
    if(CompareValue!=0)results.Branch("MinimumP2Value","TGraph",&tr,32000,0);

    results.Branch("P1Value",&p1value,"P1Value/F");
    if(CompareValue!=0)results.Branch("P2Value",&p2value,"P2Value/F");


    // Initializing
    gr->Delete();
    tr->Delete();
    gr = new TGraph();
    tr = new TGraph();

    int ia[6];
    for(ia[0]=0; ia[0]<npoints[0]; ia[0]++)  // dm32
    for(ia[1]=0; ia[1]<npoints[1]; ia[1]++)  // ssqth23
    for(ia[2]=0; ia[2]<npoints[2]; ia[2]++)  // dm21
    for(ia[3]=0; ia[3]<npoints[3]; ia[3]++)  // ssqth13
    for(ia[4]=0; ia[4]<npoints[4]; ia[4]++)  // ssqth12
    for(ia[5]=0; ia[5]<npoints[5]; ia[5]++)  // dcp
    { 

        //Clear the fitter 
        minuit->Clear(); 

	 //fParam1 = startPoints[0][0];
	 //fParam2 = startPoints[1][0];

        p1value = StandardValue;
        if(CompareValue!=0)
        {
            p2value = CompareValue;
        }
   

        int piter = 0;
  
        //Set Systematic Parameters
        for(int i=0; i<systParams->GetNSysts(); i++)
        {
            minuit->SetParameter(piter,Form("syst%d",i),systParams->GetParamNominal(i),systParams->GetParamError(i)/10.,0,0);
            if(DensityKD == 2.6 && systParams->GetParamType(i)==MATTERDENSITY)  // In case of 2 HK calculation
            {
                minuit->FixParameter(piter);	
            }
            piter++;
        }

        //Set Oscillation parameters
        for(int i=0; i<oscParams->GetNOsc(); i++)
        {
            double param_value = oscParams->GetParamValue(i);
            if(i == 1)  // Modify the sin theta 23 value
            {
                param_value = StandardValue;    
            } 
            if(oscParams->GetParamType(i)==0)   //  if the parameter has a flat prior
            {
                minuit->SetParameter(piter,Form("osc%d",i),param_value,oscParams->GetParamError(i)/100.,
                                    (i==0 && startPoints[i][ia[i]]<0. && oscParams->GetParamMinimum(i)>0. ? -1.0*oscParams->GetParamMaximum(i) : oscParams->GetParamMinimum(i)),
                                    (i==0 && startPoints[i][ia[i]]<0. && oscParams->GetParamMaximum(i)>0. ? -1.0*oscParams->GetParamMinimum(i) : oscParams->GetParamMaximum(i)));
            }
            else    //  if it has a Gaussian prior (external constraints)
            {
                minuit->SetParameter(piter,Form("osc%d",i),param_value,oscParams->GetParamError(i)/100.,oscParams->GetParamMinimum(i),oscParams->GetParamMaximum(i));
            }
            if(fixedParam[i])   // In case of fixing parameter
            {
                minuit->FixParameter(piter); // Fixing oscillation parameters we choose
            }
            piter++;
        }  
        //Set the FCN
        minuit->SetFCN(myFcn);
           
        // set print level
        double arglist[100];
        arglist[0] = 2;
        minuit->ExecuteCommand("SET PRINT",arglist,2);

        //Set the minimizer/
        arglist[0] = 50000; // number of function calls
        arglist[1] = 0.001; // tolerance
        minuit->ExecuteCommand("MIGRAD",arglist,2);
        //minuit->ExecuteCommand("MINI",arglist,2);

        //Get the fit results
        double amin,edm,errdeff;
        int nvpar, nparx;
        int conv = minuit->GetStats(amin,edm,errdeff,nvpar,nparx);

        std::cout << "Converged : " << conv << std::endl;
        std::cout << "Minimum   : " << amin << std::endl;
        //If the fit converged and the minimum is less than other starting points, save it
        if(amin<minimum && conv > 0) minimum = amin;
 
        //Fill the ttree
        //results.Fill();
             
    }

    double tempP1value, tempP2value, tempMinimum;
    //value = fSet1 * 40 + fSet2;
    if(CompareValue!=0)
    {
        CurrentRunNumber = fSet1 * 100 + fSet2;
    } 
    else 
    {
        CurrentRunNumber = fSet1;
    }

    gr->SetPoint(CurrentRunNumber, p1value, minimum);
    if(CompareValue!=0)
    {
        tr->SetPoint(CurrentRunNumber, p2value, minimum);
    }
    
    gr->GetPoint(CurrentRunNumber, tempP1value, tempMinimum);
    if(CompareValue!=0)
    {
        tr->GetPoint(CurrentRunNumber, tempP2value, tempMinimum);
    }

    if(CompareValue!=0)
    {
        std::cout << "value : " << CurrentRunNumber << " StandardValue : " << StandardValue << " CompareValue : " << CompareValue << " minimum : " << minimum << std::endl;
        std::cout << "Temp value : " << CurrentRunNumber << " StandardValue : " << tempP1value << " CompareValue : " << tempP2value << " minimum : " << tempMinimum << std::endl;
    }
    else
    {
        std::cout << "value : " << CurrentRunNumber << " StandardValue : " << StandardValue <<  " minimum : " << minimum << std::endl;
        std::cout << "Temp value : " << CurrentRunNumber << " StandardValue : " << tempP1value << " minimum : " << tempMinimum << std::endl;
    }
    //Fill the ttree
    results.Fill();

    //Write the ttree to the output file
    results.Write();
    fout->Close();
}

//Function to parse the input arguments
void ParseArgs(int argc, char **argv)
{

    std::cout << "parse args" << std::endl;
    int nargs = 1; 
    if(argc<(nargs*2+1))
    { 
        exit(1); 
    }
    for(int i = 1; i < argc; i++)
    {
        if(std::string(argv[i]) == "-c1") fCardFileHK = argv[++i];
        if(std::string(argv[i]) == "-c2") fCardFileKD = argv[++i];
        if(std::string(argv[i]) == "-co") fCardOsc = argv[++i];
        if(std::string(argv[i]) == "-cs") fCardSyst = argv[++i];
        if(std::string(argv[i]) == "-o") fOutFile = argv[++i];
        if(std::string(argv[i]) == "-fh") fFixHierarchy = atoi(argv[++i]);
        if(std::string(argv[i]) == "-s1") fSet1 = atoi(argv[++i]);
        if(std::string(argv[i]) == "-s2") fSet2 = atoi(argv[++i]);
    } 
}

//Function to calculate the likelihood
double CalcLogLikelihood(TH1D **data1Re, TH1D **data1Rmu, TH1D **pred1Re, TH1D **pred1Rmu, bool cut)
{

    double ll = 0.;

    for(int i=0; i<2; i++)
    {
        for(int j=1; j<data1Re[i]->GetNbinsX(); j++)
        {
            double erec = data1Re[i]->GetXaxis()->GetBinCenter(j);
            if(erec>1.2 && cut) 
            {
                continue;
            }
            if(erec<0.1)
            {   
                continue;
            }

            if(pred1Re[i]->GetBinContent(j)>1.0e-9 && data1Re[i]->GetBinContent(j)>1.0e-9)
            {
                ll += 2.0*(pred1Re[i]->GetBinContent(j)-data1Re[i]->GetBinContent(j)+
                  data1Re[i]->GetBinContent(j)*TMath::Log( data1Re[i]->GetBinContent(j)/pred1Re[i]->GetBinContent(j)) );
            }
            else if(data1Re[i]->GetBinContent(j)<=1.0e-9)
            {
                ll += 2.0*pred1Re[i]->GetBinContent(j);
            }
        } 
    }
    for(int i=0; i<2; i++)
    {
        for(int j=1; j<data1Rmu[i]->GetNbinsX(); j++)
        {
            double erec = data1Rmu[i]->GetXaxis()->GetBinCenter(j);
            if(erec<0.2) 
            {
                continue;
            }
            if(pred1Rmu[i]->GetBinContent(j)>1.0e-9 && data1Rmu[i]->GetBinContent(j)>1.0e-9)
            {
                ll += 2.0*(pred1Rmu[i]->GetBinContent(j)-data1Rmu[i]->GetBinContent(j)+
                  data1Rmu[i]->GetBinContent(j)*TMath::Log( data1Rmu[i]->GetBinContent(j)/pred1Rmu[i]->GetBinContent(j)) );
            }
            else if(data1Rmu[i]->GetBinContent(j)<=1.0e-9)
            {
                ll += 2.0*pred1Rmu[i]->GetBinContent(j);
            }
        }
    } 

    return ll;
}

//FCN function for the minimizer
void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{


    double matterDensity[2];
    matterDensity[0] = histsHK->GetDensity();
    if(histsKD!=NULL)
    {
        matterDensity[1] = histsKD->GetDensity();  
    } 

    //Get the updated systematic parameter values
    for(int i=0; i<systParams->GetNSysts(); i++)
    {
        systParams->SetParamValue(i,p[i]);
        if(systParams->GetParamType(i)==MATTERDENSITY)
        {
            if(systParams->DoDetector(1,i)) matterDensity[0] = p[i];
            if(systParams->DoDetector(2,i) && DensityKD != 2.6) matterDensity[1] = p[i];
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
    if(histsKD!=NULL)
    {
        //Calculate the secon detector prediction
        prob->SetBaseLine(histsKD->GetBaseline());
        prob->SetDensity(matterDensity[1]);
        histsKD->ApplyOscillations();
        histsKD->GetPredictions(pred1ReKD, pred1RmuKD,systParams);
    }

    //Get the data chi2
    double chisq = 0.;
    chisq += CalcLogLikelihood(data1ReHK,data1RmuHK,pred1ReHK, pred1RmuHK, (histsHK->GetBaseline()<300.) ); 
    if(histsKD!=NULL) 
    {
        chisq += CalcLogLikelihood(data1ReKD,data1RmuKD,pred1ReKD, pred1RmuKD, (histsKD->GetBaseline()<300.) );
    }

    //Get the oscillation parameter penalty term
    for(int i=0; i<oscParams->GetNOsc(); i++)
    {
        if(oscParams->GetParamType(i)==1)
        {
            chisq += pow( (p[i+nSyst]-oscParams->GetParamValue(i))/oscParams->GetParamError(i), 2);
        }
    }
    //Get the systematic parameter penalty term
    chisq += systParams->GetChi2();

    fval = chisq;

}

