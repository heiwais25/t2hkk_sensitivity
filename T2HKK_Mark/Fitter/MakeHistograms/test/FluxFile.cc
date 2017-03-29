#define FluxFile_cxx
#include "FluxFile.h"
#include "TSystem.h"
#include <fstream>
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include "math.h"
#include "stdlib.h"
#include "TH1D.h"
#include "OscParams.h"
#include "SystParams.h"

#include <TRandom1.h>
using namespace std;

/*=====================================================================================================================================
	
	In this file, we will treat root file in details

=====================================================================================================================================*/

FluxFile::FluxFile(char * cardFile, ProbWrapper *prob, int detid)
{
	oscProb = prob;
	cout << cardFile << endl;
	detID = (int)(pow(2,detid)+0.00001);

	ifstream infile(cardFile);
	string line;
	getline(infile, line);
	cout << line.c_str() << endl;

	// fflux = new TFile(line.c_str(),"READ");
	fflux = new TFile("../../inputs/flux_weights_1p5_1100km.root","READ");
	// fflux = new TFile();
	std::cout << "test" << std::endl;
	fluxRatios[0][0] = (TH1D*)fflux->Get("weight_numode_numu"); 
	fluxRatios[0][1] = (TH1D*)fflux->Get("weight_numode_numub"); 
	fluxRatios[0][2] = (TH1D*)fflux->Get("weight_numode_nue"); 
	fluxRatios[0][3] = (TH1D*)fflux->Get("weight_numode_nueb"); 
	fluxRatios[1][0] = (TH1D*)fflux->Get("weight_anumode_numu");
	fluxRatios[1][1] = (TH1D*)fflux->Get("weight_anumode_numub");
	fluxRatios[1][2] = (TH1D*)fflux->Get("weight_anumode_nue"); 
	fluxRatios[1][3] = (TH1D*)fflux->Get("weight_anumode_nueb");

	  cout << "Load Raw histograms" << std::endl;

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<6; j++)
		{
			getline(infile, line);
			mcWeights[i][j] = atof(line.c_str());
			std::cout << "Weights " << i << " " << j << " " << mcWeights[i][j] << std::endl;
			getline(infile, line);
			fraw[i][j] = new TFile(line.c_str());
			//fraw[i][j].OpenFile(line.c_str(),"READ");
			raw1Re[i][j][0] = (TH2D*)fraw[i][j]->Get("enu_erec_1Re_CCQE");
			raw1Re[i][j][1] = (TH2D*)fraw[i][j]->Get("enu_erec_1Re_CCnQE");
			raw1Re[i][j][2] = (TH2D*)fraw[i][j]->Get("enu_erec_1Re_NC");
			raw1Rmu[i][j][0] = (TH2D*)fraw[i][j]->Get("enu_erec_1Rmu_CCQE");
			raw1Rmu[i][j][1] = (TH2D*)fraw[i][j]->Get("enu_erec_1Rmu_CCnQE");
			raw1Rmu[i][j][2] = (TH2D*)fraw[i][j]->Get("enu_erec_1Rmu_NC");
			raw1ReOsc[i][j][0] = (TH2D*)raw1Re[i][j][0]->Clone(Form("raw1ReOsc_%d%d0",i,j));
			raw1ReOsc[i][j][1] = (TH2D*)raw1Re[i][j][1]->Clone(Form("raw1ReOsc_%d%d1",i,j));
			raw1ReOsc[i][j][2] = (TH2D*)raw1Re[i][j][2]->Clone(Form("raw1ReOsc_%d%d2",i,j));
			raw1RmuOsc[i][j][0] = (TH2D*)raw1Rmu[i][j][0]->Clone(Form("raw1RmuOsc_%d%d0",i,j));
			raw1RmuOsc[i][j][1] = (TH2D*)raw1Rmu[i][j][1]->Clone(Form("raw1RmuOsc_%d%d1",i,j));
			raw1RmuOsc[i][j][2] = (TH2D*)raw1Rmu[i][j][2]->Clone(Form("raw1RmuOsc_%d%d2",i,j));
		}
	}  

	char pol[2][10] = {"FHC","RHC"};
	char flav[6][20] = {"numu","numub","nue","nueb","numuxnue","numubxnueb"}; 
	char imode[3][10] = {"CCQE","CCnQE","NC"};

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<6; j++)
		{
			for(int k=0; k<3; k++)
			{
				TAxis *blah = raw1Re[i][j][k]->GetYaxis();
				pred1Re[i][j][k] = new TH1D(Form("pred1Re_%s_%s_%s",pol[i],flav[j],imode[k]),"",blah->GetNbins(),blah->GetXbins()->GetArray());
				blah = raw1Rmu[i][j][k]->GetYaxis();
				pred1Rmu[i][j][k] = new TH1D(Form("pred1Rmu_%s_%s_%s",pol[i],flav[j],imode[k]),"",blah->GetNbins(),blah->GetXbins()->GetArray());
			}
		}
	}  


	getline(infile, line);
	baseline = atof(line.c_str()); 
	getline(infile, line);
	density = atof(line.c_str()); 

	getline(infile, line);
	massWeight = atof(line.c_str())/22.5;
	getline(infile, line);
	potWeight[0] = 1.04*atof(line.c_str());
	getline(infile, line);
	potWeight[1] = 1.04*atof(line.c_str());

}

FluxFile::~FluxFile()
{
	
}

void FluxFile::ApplyFluxWeights()
{

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<6; j++)
		{
			int flavId = j%4;
			for(int k=0; k<3; k++)
			{

				for(int xi=1; xi<=raw1Re[i][j][k]->GetNbinsX()+1; xi++)
				{
					double enu = raw1Re[i][j][k]->GetXaxis()->GetBinCenter(xi);
					int fbin = fluxRatios[i][flavId]->FindBin(enu);
					double fweight = fluxRatios[i][flavId]->GetBinContent(fbin);
					for(int yi=1; yi<=raw1Re[i][j][k]->GetNbinsY()+1; yi++)
					{
						raw1Re[i][j][k]->SetBinContent(xi,yi,raw1Re[i][j][k]->GetBinContent(xi,yi)*fweight*mcWeights[i][j]*massWeight*potWeight[i]);
					}
				}
			

				for(int xi=1; xi<=raw1Rmu[i][j][k]->GetNbinsX()+1; xi++)
				{
					double enu = raw1Rmu[i][j][k]->GetXaxis()->GetBinCenter(xi);
					int fbin = fluxRatios[i][flavId]->FindBin(enu);
					double fweight = fluxRatios[i][flavId]->GetBinContent(fbin);
					for(int yi=1; yi<=raw1Rmu[i][j][k]->GetNbinsY()+1; yi++)
					{
						raw1Rmu[i][j][k]->SetBinContent(xi,yi,raw1Rmu[i][j][k]->GetBinContent(xi,yi)*fweight*mcWeights[i][j]*massWeight*potWeight[i]);
					}
				}
			}
		}
	}
}

void FluxFile::BuildHistograms(SystParams *systParams)
{

	double *weightmm = new double[raw1Re[0][0][0]->GetNbinsX()]; 
	double *weightmbmb = new double[raw1Re[0][0][0]->GetNbinsX()]; 
	double *weightee = new double[raw1Re[0][0][0]->GetNbinsX()]; 
	double *weightebeb = new double[raw1Re[0][0][0]->GetNbinsX()]; 
	double *weightme = new double[raw1Re[0][0][0]->GetNbinsX()]; 
	double *weightmbeb = new double[raw1Re[0][0][0]->GetNbinsX()]; 
	for(int xi=1; xi<=raw1Re[0][0][0]->GetNbinsX(); xi++)
	{
		double enu = raw1Re[0][0][0]->GetXaxis()->GetBinCenter(xi);
		weightmm[xi-1] = oscProb->GetProbNuMuNuMu(enu);
		weightmbmb[xi-1] = oscProb->GetProbNuMuBarNuMuBar(enu);
		weightee[xi-1] = oscProb->GetProbNuENuE(enu);
		weightebeb[xi-1] = oscProb->GetProbNuEBarNuEBar(enu);
		weightme[xi-1] = oscProb->GetProbNuMuNuE(enu);
		weightmbeb[xi-1] = oscProb->GetProbNuMuBarNuEBar(enu);
		//if(enu<0.7 && enu>0.5) std::cout << weightme[xi-1] << std::endl;
	}

	for(int i=0; i<2; i++)
	{
	     int pmap = pow(2,i);
	     for(int j=0; j<6; j++)
	     {
			int fmap = pow(2,j);
			for(int k=0; k<3; k++)
			{
				int xmap = pow(2,k);
				double sysweight = 1.0;

				for(int si=0; si<systParams->GetNSysts(); si++)
				{
					if( systParams->GetParamType(si)!=ERECNORM && systParams->GetParamType(si)!=ENUNORM ) continue;
					if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
					if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
					if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
					if( (detID & systParams->GetDetectorID(si))==0 ) continue;
					sysweight *= systParams->GetParamValue(i);
				}  
         //std::cout << i << " " << j << " " << k << " " << sysweight << std::endl;
 
				double *weight1 = new double[raw1Re[i][j][k]->GetNbinsX()];
				for(int xi=1; xi<=raw1Re[i][j][k]->GetNbinsX(); xi++)
				{
					weight1[xi-1] = 1.0;
					//Osc probabiltiy weights
					if(k<2)
					{
						if(j==0) weight1[xi-1] = weightmm[xi-1];
						else if(j==1) weight1[xi-1] = weightmbmb[xi-1];
						else if(j==2) weight1[xi-1] = weightee[xi-1];
						else if(j==3) weight1[xi-1] = weightebeb[xi-1];
						else if(j==4) weight1[xi-1] = weightme[xi-1];
						else if(j==5) weight1[xi-1] = weightmbeb[xi-1];
					}
					if(k==2 && j>3) weight1[xi-1] = 0.;
					weight1[xi-1] *= sysweight;
				}
				for(int yi=1; yi<=raw1Re[i][j][k]->GetNbinsY()+1; yi++)
				{
					double count = 0.;
					for(int xi=1; xi<=raw1Re[i][j][k]->GetNbinsX(); xi++) count += raw1Re[i][j][k]->GetBinContent(xi,yi)*weight1[xi-1];
					pred1Re[i][j][k]->SetBinContent(yi,count);
				}
				delete [] weight1;

				double *weight2 = new double[raw1Rmu[i][j][k]->GetNbinsX()];
				for(int xi=1; xi<=raw1Rmu[i][j][k]->GetNbinsX(); xi++)
				{
					weight2[xi-1] = 1.0;
					//Osc probabiltiy weights
					if(k<2)
					{
						if(j==0) weight2[xi-1] = weightmm[xi-1];
						else if(j==1) weight2[xi-1] = weightmbmb[xi-1];
						else if(j==2) weight2[xi-1] = weightee[xi-1];
						else if(j==3) weight2[xi-1] = weightebeb[xi-1];
						else if(j==4) weight2[xi-1] = weightme[xi-1];
						else if(j==5) weight2[xi-1] = weightmbeb[xi-1];
			           }
					if(k==2 && j>3) weight2[xi-1] = 0.;
					weight2[xi-1] *= sysweight;
				}
				for(int yi=1; yi<=raw1Rmu[i][j][k]->GetNbinsY()+1; yi++)
				{
					double count = 0.;
					for(int xi=1; xi<=raw1Rmu[i][j][k]->GetNbinsX(); xi++) count += raw1Rmu[i][j][k]->GetBinContent(xi,yi)*weight2[xi-1];
					pred1Rmu[i][j][k]->SetBinContent(yi,count);
				}
				delete [] weight2;    
			}   
		}
	}

	//Do the energy scale uncertainties

	delete [] weightmm;
	delete [] weightmbmb;
	delete [] weightee;
	delete [] weightebeb;
	delete [] weightme;
	delete [] weightmbeb;
}

void FluxFile::BuildHistograms()
{

   double *weightmm = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightmbmb = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightee = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightebeb = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightme = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightmbeb = new double[raw1Re[0][0][0]->GetNbinsX()];
   for(int xi=1; xi<=raw1Re[0][0][0]->GetNbinsX(); xi++){
     double enu = raw1Re[0][0][0]->GetXaxis()->GetBinCenter(xi);
     weightmm[xi-1] = oscProb->GetProbNuMuNuMu(enu);
     weightmbmb[xi-1] = oscProb->GetProbNuMuBarNuMuBar(enu);
     weightee[xi-1] = oscProb->GetProbNuENuE(enu);
     weightebeb[xi-1] = oscProb->GetProbNuEBarNuEBar(enu);
     weightme[xi-1] = oscProb->GetProbNuMuNuE(enu);
     weightmbeb[xi-1] = oscProb->GetProbNuMuBarNuEBar(enu);
     //if(enu<0.7 && enu>0.5) std::cout << weightme[xi-1] << std::endl;
   }
    

   for(int i=0; i<2; i++)
     for(int j=0; j<6; j++)
       for(int k=0; k<3; k++){
         double *weight1 = new double[raw1Re[i][j][k]->GetNbinsX()];
         for(int xi=1; xi<=raw1Re[i][j][k]->GetNbinsX(); xi++){
           //double enu = raw1Re[i][j][k]->GetXaxis()->GetBinCenter(xi);
           weight1[xi-1] = 1.0;
           if(k<2){
             if(j==0) weight1[xi-1] = weightmm[xi-1];
             else if(j==1) weight1[xi-1] = weightmbmb[xi-1];
             else if(j==2) weight1[xi-1] = weightee[xi-1];
             else if(j==3) weight1[xi-1] = weightebeb[xi-1];
             else if(j==4) weight1[xi-1] = weightme[xi-1];
             else if(j==5) weight1[xi-1] = weightmbeb[xi-1];
           }  
           if(k==2 && j>3) weight1[xi-1] = 0.;
         }  
         for(int yi=1; yi<=raw1Re[i][j][k]->GetNbinsY()+1; yi++){
           double count = 0.;
           for(int xi=1; xi<=raw1Re[i][j][k]->GetNbinsX(); xi++) count += raw1Re[i][j][k]->GetBinContent(xi,yi)*weight1[xi-1];
           pred1Re[i][j][k]->SetBinContent(yi,count);
         }
         delete [] weight1;

         double *weight2 = new double[raw1Rmu[i][j][k]->GetNbinsX()];
         for(int xi=1; xi<=raw1Rmu[i][j][k]->GetNbinsX(); xi++){
           //double enu = raw1Rmu[i][j][k]->GetXaxis()->GetBinCenter(xi);
           weight2[xi-1] = 1.0;
           if(k<2){
             if(j==0) weight2[xi-1] = weightmm[xi-1];
             else if(j==1) weight2[xi-1] = weightmbmb[xi-1];
             else if(j==2) weight2[xi-1] = weightee[xi-1];
             else if(j==3) weight2[xi-1] = weightebeb[xi-1];
             else if(j==4) weight2[xi-1] = weightme[xi-1];
             else if(j==5) weight2[xi-1] = weightmbeb[xi-1];
           }  
           if(k==2 && j>3) weight2[xi-1] = 0.;
         }  
         for(int yi=1; yi<=raw1Rmu[i][j][k]->GetNbinsY()+1; yi++){
           double count = 0.;
           for(int xi=1; xi<=raw1Rmu[i][j][k]->GetNbinsX(); xi++) count += raw1Rmu[i][j][k]->GetBinContent(xi,yi)*weight2[xi-1];
           pred1Rmu[i][j][k]->SetBinContent(yi,count);
         } 
         delete [] weight2;   
      }

   delete [] weightmm;
   delete [] weightmbmb;
   delete [] weightee;
   delete [] weightebeb;
   delete [] weightme;
   delete [] weightmbeb;

}



void FluxFile::ApplyOscillations()
{

   double *weightmm = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightmbmb = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightee = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightebeb = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightme = new double[raw1Re[0][0][0]->GetNbinsX()];
   double *weightmbeb = new double[raw1Re[0][0][0]->GetNbinsX()];
   for(int xi=1; xi<=raw1Re[0][0][0]->GetNbinsX(); xi++){
     double enu = raw1Re[0][0][0]->GetXaxis()->GetBinCenter(xi);
     weightmm[xi-1] = oscProb->GetProbNuMuNuMu(enu);
     weightmbmb[xi-1] = oscProb->GetProbNuMuBarNuMuBar(enu);
     weightee[xi-1] = oscProb->GetProbNuENuE(enu);
     weightebeb[xi-1] = oscProb->GetProbNuEBarNuEBar(enu);
     weightme[xi-1] = oscProb->GetProbNuMuNuE(enu);
     weightmbeb[xi-1] = oscProb->GetProbNuMuBarNuEBar(enu);
     //if(enu<0.7 && enu>0.5) std::cout << weightme[xi-1] << std::endl;
   }
    

   for(int i=0; i<2; i++)
     for(int j=0; j<6; j++){
       double *weight1;
       if(j==0) weight1 = weightmm;
       else if(j==1) weight1 = weightmbmb;
       else if(j==2) weight1 = weightee;
       else if(j==3) weight1 = weightebeb;
       else if(j==4) weight1 = weightme;
       else if(j==5) weight1 = weightmbeb;
       for(int k=0; k<3; k++){

         for(int xi=1; xi<=raw1Re[i][j][k]->GetNbinsX(); xi++){
           double w = weight1[xi-1];
           if(k==2) w = 1.0;
           if(k==2 && j>3) w = 0.0;
           for(int yi=1; yi<=raw1Re[i][j][k]->GetNbinsY(); yi++) raw1ReOsc[i][j][k]->SetBinContent(xi,yi,w*raw1Re[i][j][k]->GetBinContent(xi,yi));
         }   

         for(int xi=1; xi<=raw1Rmu[i][j][k]->GetNbinsX(); xi++){
           double w = weight1[xi-1];
           if(k==2) w = 1.0;
           if(k==2 && j>3) w = 0.0;
           for(int yi=1; yi<=raw1Rmu[i][j][k]->GetNbinsY(); yi++) raw1RmuOsc[i][j][k]->SetBinContent(xi,yi,w*raw1Rmu[i][j][k]->GetBinContent(xi,yi));
         }   

       }
     }  

   delete [] weightmm;
   delete [] weightmbmb;
   delete [] weightee;
   delete [] weightebeb;
   delete [] weightme;
   delete [] weightmbeb;

}

void FluxFile::SaveToFile(char *filename)
{
 
   TFile *fout = new TFile(filename,"RECREATE");
   for(int i=0; i<2; i++)
     for(int j=0; j<6; j++)
       for(int k=0; k<3; k++){
         pred1Re[i][j][k]->Write();
         pred1Rmu[i][j][k]->Write();
       }
   fout->Close();

}  

void FluxFile::GetPredictions(TH1D **hists1Re, TH1D **hists1Rmu)
{

  /*int nbins = 56;
  double bins[57] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,
                     1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
                     3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4,4.6,4.8,
                     5.0,5.2,5.4,5.6,5.8,6.0,6.5,7.0,7.5,8.0,9.0,10.0};*/

  if(hists1Re[0]!=NULL) hists1Re[0]->Delete();
  if(hists1Re[1]!=NULL) hists1Re[1]->Delete();
  if(hists1Rmu[0]!=NULL) hists1Rmu[0]->Delete();
  if(hists1Rmu[1]!=NULL) hists1Rmu[1]->Delete();

  hists1Re[0] = (TH1D*)pred1Re[0][0][0]->Clone("FHC1Re");
  hists1Re[1] = (TH1D*)pred1Re[1][0][0]->Clone("RHC1Re");
  hists1Rmu[0] = (TH1D*)pred1Rmu[0][0][0]->Clone("FHC1Rmu");
  hists1Rmu[1] = (TH1D*)pred1Rmu[1][0][0]->Clone("RHC1Rmu");

  /*hists1Re[0] = (TH1D*)pred1Re[0][0][0]->Rebin(nbins,"FHC1Re",bins);
  hists1Re[1] = (TH1D*)pred1Re[1][0][0]->Rebin(nbins,"RHC1Re",bins);
  hists1Rmu[0] = (TH1D*)pred1Rmu[0][0][0]->Rebin(nbins,"FHC1Rmu",bins);
  hists1Rmu[1] = (TH1D*)pred1Rmu[1][0][0]->Rebin(nbins,"RHC1Rmu",bins);*/

  for(int i=0; i<2; i++)
     for(int j=0; j<6; j++)
       for(int k=0; k<3; k++){
         if(j==0 && k==0) continue;
         hists1Re[i]->Add(pred1Re[i][j][k],1.0);
         hists1Rmu[i]->Add(pred1Rmu[i][j][k],1.0);
       }

}


void FluxFile::GetPredictions(TH1D **hists1Re, TH1D **hists1Rmu, SystParams *systParams, double energyScale)
{

  if(hists1Re[0]!=NULL) hists1Re[0]->Delete();
  if(hists1Re[1]!=NULL) hists1Re[1]->Delete();
  if(hists1Rmu[0]!=NULL) hists1Rmu[0]->Delete();
  if(hists1Rmu[1]!=NULL) hists1Rmu[1]->Delete();

  hists1Re[0] = (TH1D*)pred1Re[0][0][0]->Clone("FHC1Re");
  hists1Re[1] = (TH1D*)pred1Re[1][0][0]->Clone("RHC1Re");
  hists1Rmu[0] = (TH1D*)pred1Rmu[0][0][0]->Clone("FHC1Rmu");
  hists1Rmu[1] = (TH1D*)pred1Rmu[1][0][0]->Clone("RHC1Rmu");

  //Make histograms with systematics included
  for(int i=0; i<2; i++){
     int pmap = pow(2,i);
     for(int j=0; j<6; j++){
       int fmap = pow(2,j);
       for(int k=0; k<3; k++){
         int xmap = pow(2,k);
         //double sysweight = 1.0;

         /*for(int si=0; si<systParams->GetNSysts(); si++){
           if( systParams->GetParamType(si)!=ERECNORM && systParams->GetParamType(si)!=ENUNORM ) continue;
           if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
           if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
           if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
           if( (detID & systParams->GetDetectorID(si))==0 ) continue;
           sysweight *= systParams->GetParamValue(i);
         }*/

         double *weightErec = new double[raw1ReOsc[i][j][k]->GetNbinsY()];
         double *weightEnu = new double[raw1ReOsc[i][j][k]->GetNbinsX()];
         int pidmap = 1;

         //Save the Erec weights
         for(int yi=1; yi<=raw1ReOsc[i][j][k]->GetNbinsY(); yi++){
           weightErec[yi-1] = 1.0;
           double erec = raw1ReOsc[i][j][k]->GetYaxis()->GetBinCenter(yi);
           for(int si=0; si<systParams->GetNSysts(); si++){
             if( systParams->GetParamType(si)!=ERECNORM) continue;
             if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
             if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
             if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
             if( (detID & systParams->GetDetectorID(si))==0 ) continue;
             if( (pidmap & systParams->GetPIDType(si))==0 ) continue;
             if( erec<systParams->GetMaxEnergy(si) && erec>systParams->GetMinEnergy(si) ) weightErec[yi-1] *= systParams->GetParamValue(si);
           }            
         }  

         //Save the Enu weights
         for(int xi=1; xi<=raw1ReOsc[i][j][k]->GetNbinsX(); xi++){
           weightEnu[xi-1] = 1.0;
           double enu = raw1ReOsc[i][j][k]->GetXaxis()->GetBinCenter(xi);
           for(int si=0; si<systParams->GetNSysts(); si++){
             if( systParams->GetParamType(si)!=ENUNORM) continue;
             if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
             if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
             if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
             if( (detID & systParams->GetDetectorID(si))==0 ) continue;
             if( (pidmap & systParams->GetPIDType(si))==0 ) continue;
             if( enu<systParams->GetMaxEnergy(si) && enu>systParams->GetMinEnergy(si) ) weightEnu[xi-1] *= systParams->GetParamValue(si);
           }
         }


         for(int yi=1; yi<=raw1ReOsc[i][j][k]->GetNbinsY(); yi++){
           double count = 0.;
           for(int xi=1; xi<=raw1ReOsc[i][j][k]->GetNbinsX(); xi++) count += raw1ReOsc[i][j][k]->GetBinContent(xi,yi)*weightErec[yi-1]*weightEnu[xi-1];
           if(j==0 && k==0) hists1Re[i]->SetBinContent(yi,count);
           else hists1Re[i]->SetBinContent(yi,count+hists1Re[i]->GetBinContent(yi));
         }

         delete [] weightErec;
         delete [] weightEnu;
         weightErec = new double[raw1RmuOsc[i][j][k]->GetNbinsY()];
         weightEnu = new double[raw1RmuOsc[i][j][k]->GetNbinsX()];
         pidmap = 2;

         //Save the Erec weights
         for(int yi=1; yi<=raw1RmuOsc[i][j][k]->GetNbinsY(); yi++){
           weightErec[yi-1] = 1.0;
           double erec = raw1RmuOsc[i][j][k]->GetYaxis()->GetBinCenter(yi);
           for(int si=0; si<systParams->GetNSysts(); si++){
             if( systParams->GetParamType(si)!=ERECNORM) continue;
             if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
             if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
             if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
             if( (detID & systParams->GetDetectorID(si))==0 ) continue;
             if( (pidmap & systParams->GetPIDType(si))==0 ) continue;
             if( erec<systParams->GetMaxEnergy(si) && erec>systParams->GetMinEnergy(si) ) weightErec[yi-1] *= systParams->GetParamValue(si);
           }
         }

         //Save the Enu weights
         for(int xi=1; xi<=raw1RmuOsc[i][j][k]->GetNbinsX(); xi++){
           weightEnu[xi-1] = 1.0;
           double enu = raw1RmuOsc[i][j][k]->GetXaxis()->GetBinCenter(xi);
           for(int si=0; si<systParams->GetNSysts(); si++){
             if( systParams->GetParamType(si)!=ENUNORM) continue;
             if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
             if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
             if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
             if( (detID & systParams->GetDetectorID(si))==0 ) continue;
             if( (pidmap & systParams->GetPIDType(si))==0 ) continue;
             if( enu<systParams->GetMaxEnergy(si) && enu>systParams->GetMinEnergy(si) ) weightEnu[xi-1] *= systParams->GetParamValue(si);
           }
         }
 
          
         for(int yi=1; yi<=raw1RmuOsc[i][j][k]->GetNbinsY(); yi++){
           double count = 0.;
           for(int xi=1; xi<=raw1RmuOsc[i][j][k]->GetNbinsX(); xi++) count += raw1RmuOsc[i][j][k]->GetBinContent(xi,yi)*weightErec[yi-1]*weightEnu[xi-1];
           if(j==0 && k==0) hists1Rmu[i]->SetBinContent(yi,count);
           else hists1Rmu[i]->SetBinContent(yi,count+hists1Rmu[i]->GetBinContent(yi));
         }

         delete [] weightErec;
         delete [] weightEnu;
       }
     }
   }         

 //Do the energy scale systematics
  for(int i=0; i<2; i++){
   int pmap = pow(2,i);
   for(int si=0; si<systParams->GetNSysts(); si++){
     if( systParams->GetParamType(si)!=ESCALE ) continue;
     if( (detID & systParams->GetDetectorID(si))==0 ) continue;
     if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
     double escale = systParams->GetParamValue(si);
     escale = energyScale;
     //Make a temporary Graph()
     double tmp = 0.;

     int pidmap = 1;
     double gre_xvalue;
     double gre_yvalue;
     if( (pidmap & systParams->GetPIDType(si))!=0 ){
       TGraph *gre = new TGraph();
       for(int l=1; l<=hists1Re[i]->GetNbinsX(); l++)
       {
         gre->SetPoint(l-1,hists1Re[i]->GetXaxis()->GetBinCenter(l)*escale,hists1Re[i]->GetBinContent(l)/(hists1Re[i]->GetXaxis()->GetBinWidth(l)*escale));
         // escale_nue->SetPoint(l-1,hists1Re[i]->GetXaxis()->GetBinCenter(l)*escale,hists1Re[i]->GetBinContent(l)/(hists1Re[i]->GetXaxis()->GetBinWidth(l)*escale));
       }

       for(int l=1; l<=hists1Re[i]->GetNbinsX(); l++){
         tmp = gre->Eval(hists1Re[i]->GetXaxis()->GetBinCenter(l),0,"S");
         // std::cout << "tmp value" << tmp << std::endl;
         if(tmp<0.) tmp = 0.;
         hists1Re[i]->SetBinContent(l,tmp*hists1Re[i]->GetXaxis()->GetBinWidth(l));
       } 
       gre->Delete();
     }

     pidmap = 2;
     if( (pidmap & systParams->GetPIDType(si))!=0 ){
       TGraph *grmu = new TGraph();
       for(int l=1; l<=hists1Rmu[i]->GetNbinsX(); l++)
       {
         grmu->SetPoint(l-1,hists1Rmu[i]->GetXaxis()->GetBinCenter(l)*escale,hists1Rmu[i]->GetBinContent(l)/(hists1Rmu[i]->GetXaxis()->GetBinWidth(l)*escale));
         // escale_numu->SetPoint(l-1,hists1Rmu[i]->GetXaxis()->GetBinCenter(l)*escale,hists1Rmu[i]->GetBinContent(l)/(hists1Rmu[i]->GetXaxis()->GetBinWidth(l)*escale));

       }
       for(int l=1; l<=hists1Rmu[i]->GetNbinsX(); l++){
         tmp = grmu->Eval(hists1Rmu[i]->GetXaxis()->GetBinCenter(l),0,"S");
         if(tmp<0.) tmp = 0.;
         hists1Rmu[i]->SetBinContent(l,tmp*hists1Rmu[i]->GetXaxis()->GetBinWidth(l));
       } 
       grmu->Delete();
     }  
   }
 }  
 

}

void FluxFile::BuildHistogramsSystOnly(SystParams *systParams )
{

  //Make histograms with systematics included
  for(int i=0; i<2; i++){
     int pmap = pow(2,i);
     for(int j=0; j<6; j++){
       int fmap = pow(2,j);
       for(int k=0; k<3; k++){
         int xmap = pow(2,k);
         //double sysweight = 1.0;

         double *weightErec = new double[raw1Re[i][j][k]->GetNbinsY()];
         double *weightEnu = new double[raw1Re[i][j][k]->GetNbinsX()];
         int pidmap = 1;

         //Save the Erec weights
         for(int yi=1; yi<=raw1Re[i][j][k]->GetNbinsY(); yi++){
           weightErec[yi-1] = 1.0;
           double erec = raw1Re[i][j][k]->GetYaxis()->GetBinCenter(yi);
           for(int si=0; si<systParams->GetNSysts(); si++){
             if( systParams->GetParamType(si)!=ERECNORM) continue;
             if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
             if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
             if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
             if( (detID & systParams->GetDetectorID(si))==0 ) continue;
             if( (pidmap & systParams->GetPIDType(si))==0 ) continue;
             if( erec<systParams->GetMaxEnergy(si) && erec>systParams->GetMinEnergy(si) ) weightErec[yi-1] *= systParams->GetParamValue(si);
           }            
         }  

         //Save the Enu weights
         for(int xi=1; xi<=raw1Re[i][j][k]->GetNbinsX(); xi++){
           weightEnu[xi-1] = 1.0;
           double enu = raw1Re[i][j][k]->GetXaxis()->GetBinCenter(xi);
           for(int si=0; si<systParams->GetNSysts(); si++){
             if( systParams->GetParamType(si)!=ENUNORM) continue;
             if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
             if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
             if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
             if( (detID & systParams->GetDetectorID(si))==0 ) continue;
             if( (pidmap & systParams->GetPIDType(si))==0 ) continue;
             if( enu<systParams->GetMaxEnergy(si) && enu>systParams->GetMinEnergy(si) ) weightEnu[xi-1] *= systParams->GetParamValue(si);
           }
         }


         for(int yi=1; yi<=raw1Re[i][j][k]->GetNbinsY(); yi++)
           for(int xi=1; xi<=raw1Re[i][j][k]->GetNbinsX(); xi++) 
             pred1Re2D[i][j][k]->SetBinContent(xi,yi,raw1Re[i][j][k]->GetBinContent(xi,yi)*weightErec[yi-1]*weightEnu[xi-1]);
            
         delete [] weightErec;
         delete [] weightEnu;

         weightErec = new double[raw1Rmu[i][j][k]->GetNbinsY()];
         weightEnu = new double[raw1Rmu[i][j][k]->GetNbinsX()];
         pidmap = 2;

         //Save the Erec weights
         for(int yi=1; yi<=raw1Rmu[i][j][k]->GetNbinsY(); yi++){
           weightErec[yi-1] = 1.0;
           double erec = raw1Rmu[i][j][k]->GetYaxis()->GetBinCenter(yi);
           for(int si=0; si<systParams->GetNSysts(); si++){
             if( systParams->GetParamType(si)!=ERECNORM) continue;
             if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
             if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
             if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
             if( (detID & systParams->GetDetectorID(si))==0 ) continue;
             if( (pidmap & systParams->GetPIDType(si))==0 ) continue;
             if( erec<systParams->GetMaxEnergy(si) && erec>systParams->GetMinEnergy(si) ) weightErec[yi-1] *= systParams->GetParamValue(si);
           }
         }

         //Save the Enu weights
         for(int xi=1; xi<=raw1Rmu[i][j][k]->GetNbinsX(); xi++){
           weightEnu[xi-1] = 1.0;
           double enu = raw1Rmu[i][j][k]->GetXaxis()->GetBinCenter(xi);
           for(int si=0; si<systParams->GetNSysts(); si++){
             if( systParams->GetParamType(si)!=ENUNORM) continue;
             if( (fmap & systParams->GetFluxMode(si))==0 ) continue;
             if( (xmap & systParams->GetXsecMode(si))==0 ) continue;
             if( (pmap & systParams->GetPolarityMode(si))==0 ) continue;
             if( (detID & systParams->GetDetectorID(si))==0 ) continue;
             if( (pidmap & systParams->GetPIDType(si))==0 ) continue;
             if( enu<systParams->GetMaxEnergy(si) && enu>systParams->GetMinEnergy(si) ) weightEnu[xi-1] *= systParams->GetParamValue(si);
           }
         }
 
          
         for(int yi=1; yi<=raw1Rmu[i][j][k]->GetNbinsY(); yi++)
           for(int xi=1; xi<=raw1Rmu[i][j][k]->GetNbinsX(); xi++) 
             pred1Rmu2D[i][j][k]->SetBinContent(xi,yi,raw1Rmu[i][j][k]->GetBinContent(xi,yi)*weightErec[yi-1]*weightEnu[xi-1]);

         delete [] weightErec;
         delete [] weightEnu;
       }
     }
   }         

}

TH2D* FluxFile::GetWeightedTemplate(FLAVOR nuflavor, RCODE reactcode, POLARITY hornpolarity, SELECTION sample)
{

  if(sample==ONERINGMU)
    return pred1Rmu2D[hornpolarity][nuflavor][reactcode];
  else
    return pred1Re2D[hornpolarity][nuflavor][reactcode];
}




int main(void)
{
	char fCardOsc[] = "../../osc_param.card";
	char fCardSyst[] = "../../syst_param_total.card";
	gSystem->Load("libTree");
	FluxFile *histsHK=NULL;
	FluxFile *histsKD=NULL;

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
	int ithier = 0;

	prob->SetDeltaCP(0);
	prob->SetBaseLine(1100.);
	prob->SetDensity(3.0);

	//Load the oscillation parameters
	oscParams = new OscParams(fCardOsc);

	//Load the systematic parameters 
	systParams = new SystParams(fCardSyst);


	// Treating data
	char cardFile[] = "../../input_kd_1p5_test.card";
	prob->SetOscillationParameters((ithier==0 ? 1.0 : -1.0)*oscParams->GetParamValue(0),oscParams->GetParamValue(1),oscParams->GetParamValue(2),oscParams->GetParamValue(3),oscParams->GetParamValue(4),oscParams->GetParamValue(5));
	FluxFile * ffile = new FluxFile(cardFile, prob, 0);
	ffile -> ApplyFluxWeights();
	ffile -> ApplyOscillations();

	TH1D * fluxfile = NULL;
	TH2D * rawFile = NULL;

	rawFile = ffile -> Getraw1Re(0, 0, 0);
	fluxfile = ffile -> GetfluxRatios(0, 0);

	for(int xi=1; xi<=rawFile->GetNbinsX()+1;xi++)
	{
		double enu = rawFile->GetXaxis()->GetBinCenter(xi);
		int fbin = fluxfile->FindBin(enu);
		double fweight = fluxfile->GetBinContent(fbin);
		// cout << enu << "\t" << log(enu) << "\t" << fbin << "\t" << fweight << endl;
		cout << rawFile->GetXaxis()->GetBinWidth(xi) << endl;
	}

	systParams->ResetParams();
	ffile->GetPredictions(data1ReKD, data1RmuKD, systParams, 0.8);



	/*===========================================================================================================================================
	 
		Deal with data before and after applying energy scale parameter
	
	 ===========================================================================================================================================*/
	char fOutFile[] = "output0.8.root";
	TFile *fout = new TFile(fOutFile,"RECREATE");
	for(int i = 0; i < 2; i++)
	{
		data1ReKD[i]->Write();
		data1RmuKD[i]->Write();
	}
	fout->Close();


	return 0;

}