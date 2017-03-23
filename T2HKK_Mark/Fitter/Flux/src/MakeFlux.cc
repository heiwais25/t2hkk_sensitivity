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

	double * carrier;

	ifstream infile(fluxFile.c_str());
	string line;
	StringParse * strParse = new StringParse();

	int Contents = 0;
	int i = 0, length = 0;
	int zenith_Count  	= 0;
	int azimuth_Count 	= -1; // To start 0 when we use -1
	int contents_Count 	= 0;
	int signal;
	while(getline(infile, line))
	{
		
		signal 	= strParse->NumberExtract(line);
		carrier 	= strParse->NumberGet();
		if(signal == 0)
		{
			// cout << "This line is blank line" << endl;
		}
		else if(signal == 1)
		{
			azimuth_Count++;
			contents_Count = 0;

			if(azimuth_Count == 12) // Initializing
			{
				azimuth_Count = 0;
				zenith_Count++;
			}
		}	
		else if(signal == 2)
		{	
			// zenith_Count 		= 0 : 0.90 ~ 1.00 / 19 : -1.0 ~ -0.90
			// azimuth_Count 		= 0 : 0 ~ 30 / 11 : 330 ~ 360
			// contents_Count 	= from 0 to 100 (total line is 101)
			hondaEnergy[zenith_Count][azimuth_Count][contents_Count] 	= carrier[0];
			hondaNuMu[zenith_Count][azimuth_Count][contents_Count]	= carrier[1];
			hondaNuMuBar[zenith_Count][azimuth_Count][contents_Count] 	= carrier[2];
			hondaNuE[zenith_Count][azimuth_Count][contents_Count] 		= carrier[3];
			hondaNuEBar[zenith_Count][azimuth_Count][contents_Count]	= carrier[4];
			contents_Count++;			
		}
	}
	int check = 19;
	int check2 = 10;
	for(int i = 0; i < 101; i++)
	{
		cout << hondaEnergy[check][check2][i] << "\t" << hondaNuMu[check][check2][i] << "\t" << hondaNuMuBar[check][check2][i] 
		<< "\t" << hondaNuE[check][check2][i] << "\t" << hondaNuEBar[check][check2][i] << endl;

	}
}

void MakeFlux::Test()
{
}


int main(void)
{
	/* code */
	string cardFile = "../resources/honda-2015-spl-solmin.d";
	MakeFlux * flux = new MakeFlux(cardFile);
	flux->Test();
	return 0;
}