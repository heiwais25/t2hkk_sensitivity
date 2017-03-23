#ifndef MakeFlux_h
#define MakeFlux_h
#define Contents_Length 101
// #include <TROOT.h>
// #include <TChain.h>
// #include <TFile.h>
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <vector>
// #include "TH2D.h"
// #include "probWrapper.h"
// #include "SystParams.h"
#include <string>
using namespace std;
class MakeFlux
{
	public:
		MakeFlux(string fluxFile);
		void Test();

	private:
		string flavor[4];
		double hondaEnergy[20][12][101];
		double hondaNuMu[20][12][101];
		double hondaNuMuBar[20][12][101];
		double hondaNuE[20][12][101];
		double hondaNuEBar[20][12][101];
};

#endif