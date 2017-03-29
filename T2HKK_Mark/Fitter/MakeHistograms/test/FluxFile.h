#ifndef FluxFile_h
#define FluxFile_h

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
#include <TGraph.h>

// #include "SystParams.h"

enum FLAVOR {NUMU, NUMUBAR, NUE, NUEBAR, NUMUNUE, NUMUBARNUEBAR};
enum RCODE {CCQE, CCnQE, NC};
enum POLARITY {FHC, RHC};
enum SELECTION {ONERINGMU, ONERINGE};

class FluxFile
{
	public:
		FluxFile(char * cardFile, ProbWrapper *prob, int detid);
		~FluxFile();
		void ApplyFluxWeights();
		void BuildHistograms(SystParams *systParams);
		void BuildHistograms();
		void BuildHistogramsSystOnly(SystParams *systParams );
		void ApplyOscillations();
		void SaveToFile(char *filename);
		void GetPredictions(TH1D **hists1Re, TH1D **hists1Rmu);
		void GetPredictions(TH1D **hists1Re, TH1D **hists1Rmu, SystParams *systParams, double energyScale);
		TH2D* GetWeightedTemplate(FLAVOR nuflavor, RCODE reactcode, POLARITY hornpolarity, SELECTION sample);

		double GetBaseline(){return baseline;};
		double GetDensity(){return density;};

		TH2D * Getraw1Re(int polarity, int flavor, int interaction){return raw1Re[polarity][flavor][interaction];};
		TH2D * Getraw1Rmu(int polarity, int flavor, int interaction){return raw1Rmu[polarity][flavor][interaction];};
		TH1D * Getpred1Re(int polarity, int flavor, int interaction){return pred1Re[polarity][flavor][interaction];};
		TH1D * Getpred1Rmu(int polarity, int flavor, int interaction){return pred1Rmu[polarity][flavor][interaction];};
		TH2D * Getraw1ReOsc(int polarity, int flavor, int interaction){return raw1ReOsc[polarity][flavor][interaction];};
		TH2D * Getraw1RmuOsc(int polarity, int flavor, int interaction){return raw1RmuOsc[polarity][flavor][interaction];};
		TH1D * GetfluxRatios(int polarity, int flavor){return fluxRatios[polarity][flavor];};
		TGraph * GetEscaleNue(){return escale_nue;};
		TGraph * GetEscaleNumu(){return escale_numu;};

	private:
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
		TFile * fflux;
		TFile * fraw[2][6];
		double baseline;
		double density;
		TGraph * escale_nue;
		TGraph * escale_numu;
		int detID;

};


#endif