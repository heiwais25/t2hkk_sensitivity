#ifndef MakeHistograms_h
#define MakeHistograms_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TH2D.h"
#include "probWrapper.h"
#include "SystParams.h"

enum FLAVOR {NUMU, NUMUBAR, NUE, NUEBAR, NUMUNUE, NUMUBARNUEBAR};
enum RCODE {CCQE, CCnQE, NC};
enum POLARITY {FHC, RHC};
enum SELECTION {ONERINGMU, ONERINGE};

class MakeHistograms {

  public :

   MakeHistograms(char *cardFile, ProbWrapper *prob, int detid=0);
   //MakeHistograms(char *cardFile, OscProb *prob);
   ~MakeHistograms();
   void ApplyFluxWeights();
   void BuildHistograms();
   void ApplyOscillations();
   void BuildHistograms(SystParams *systParams);
   void BuildHistogramsSystOnly(SystParams *systParams );
   void SaveToFile(char *filename);
   void GetPredictions(TH1D **hists1Re, TH1D **hists1Rmu);
   void GetPredictions(TH1D **hists1Re, TH1D **hists1Rmu, SystParams *systParams);
   double GetBaseline(){return baseline;};
   double GetDensity(){return density;};
   void SetMassWeight(double weight){massWeight = weight;};
   void SetPOTWeight(double fhcweight, double rhcweight){potWeight[0]=fhcweight; potWeight[1]=rhcweight;};
   void GetPOTWeight(double * inputWeight0, double * inputWeight1){*inputWeight0 = potWeight[0]; *inputWeight1 = potWeight[1];};
   TH2D* GetWeightedTemplate(FLAVOR nuflavor, RCODE reactcode, POLARITY hornpolarity, SELECTION sample);

  private :

   TH1D *fluxRatios[2][4];
   TH1D *pred1Re[2][6][3];
   TH1D *pred1Rmu[2][6][3];
   TH2D *raw1Re[2][6][3];
   TH2D *raw1Rmu[2][6][3];
   TH2D *raw1ReOsc[2][6][3];
   TH2D *raw1RmuOsc[2][6][3];
   TH2D *pred1Re2D[2][6][3];
   TH2D *pred1Rmu2D[2][6][3];
   ProbWrapper *oscProb;
   //OscProb *oscProb;
   double mcWeights[2][6];
   double massWeight;
   double potWeight[2];
   TFile *fflux;
   TFile *fraw[2][6];
   double baseline;
   double density;
   int detID;

};

#endif

