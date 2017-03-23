// #define MakeFlux_cxx
// #include <fstream>
// #include <string>
// #include <TH2.h>
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TMath.h>
// #include <TGraph.h>
// #include "math.h"
// #include "stdlib.h"
// #include "TH1D.h"
#include "MakeFlux.h"
#include "StringParse.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <stdio.h>
#include <sstream>
using namespace std;

MakeFlux::MakeFlux(string fluxFile)
{
	/* flavor declaration */
	flavor[0] = "numu";
	flavor[1] = "numubar";
	flavor[2] = "nue";
	flavor[3] = "nuebar";

	double * Angle;

	ifstream infile(fluxFile.c_str());
	string line;
	StringParse * strParse = new StringParse();
	int Contents = 0;
	int i = 0, length = 0;
	double carrier[10];
	int zenith_Count = 0;
	int azimuth_Count = -1; // To start 0 when we use -1
	int contents_Count = 0;
	double tmpValue;

	while(getline(infile, line))
	{
		if(line.at(0) == 'a')
		{
			
			strParse->NumberExtract(line);
			Angle = strParse->NumberGet();
			cout << "Cosine Zenith angle and Azimuth angle" << endl;
			cout << Angle[0] << "\t\t" << Angle[1] << "\t\t" << Angle[2] << "\t\t" << Angle[3] << endl;
			Contents = 1;
			azimuth_Count++;
		}
		else if(Contents == 1)
		{
			length = 0;
			stringstream linestream(line);
			while(linestream >> tmpValue)
			{
				if(tmpValue != 0)
				{
					carrier[length] = tmpValue;
					cout << carrier[length] << "\t\t";
					length++;
				}
			}
			if(length != 0)
			{
				// cout << carrier[0] << endl;
				cout << endl;

				hondaEnergy[zenith_Count][azimuth_Count][contents_Count] 	= carrier[0];
				hondaNuMu[zenith_Count][azimuth_Count][contents_Count]	= carrier[1];
				hondaNuMuBar[zenith_Count][azimuth_Count][contents_Count] 	= carrier[2];
				hondaNuE[zenith_Count][azimuth_Count][contents_Count] 		= carrier[3];
				hondaNuEBar[zenith_Count][azimuth_Count][contents_Count]	= carrier[4];
				// cout << hondaEnergy[zenith_Count][azimuth_Count][contents_Count] << endl;
				contents_Count++;

			}

			if(azimuth_Count == 12) // Initializing
			{
				azimuth_Count 	= 0;
				contents_Count 	= 0;
				zenith_Count++;
			}
		}
		// i++;
		// if(i == 10000)
		// {
		// 	break;
		// }
	}

	cout << endl;
	int check = 5;
	for(int i = 0; i < 101; i++)
	{
		cout << hondaEnergy[check][0][i] << "\t" << hondaNuMu[check][0][i] << "\t" << hondaNuMuBar[check][0][i] 
		<< "\t" << hondaNuE[check][0][i] << "\t" << hondaNuEBar[check][0][i] << endl;

	}


}

void MakeFlux::Test()
{
	// cout << "Hello" << endl;
	// cout << flavor[0] << endl;
}


int main(void)
{
	/* code */
	string cardFile = "../resources/honda-2015-spl-solmin.d";
	MakeFlux * flux = new MakeFlux(cardFile);
	flux->Test();
	return 0;
}