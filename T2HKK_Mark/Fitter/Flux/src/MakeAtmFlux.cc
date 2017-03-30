#include "MakeAtmFlux.h"
#include "StringParse.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <stdio.h>
#include <sstream>
#include <math.h>
using namespace std;

MakeAtmFlux::MakeAtmFlux(string fluxFile)
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

	for(int i = 0; i < CosZenith_Length; i++)
	{
		cosineZenith[i] = 0.95 - 0.10 * (double)i;
	}
	int check = 19;
	int check2 = 10;
	for(int i = 0; i < 101; i++)
	{
		// cout << log10(hondaEnergy[check][check2][i]) << "\t" << hondaNuMu[check][check2][i] << "\t" << hondaNuMuBar[check][check2][i] 
		// << "\t" << hondaNuE[check][check2][i] << "\t" << hondaNuEBar[check][check2][i] << endl;
		// cout << log10(hondaEnergy[check][check2][i]) << "\t" << log10(hondaNuMu[check][check2][i]) << endl;

	}
	
}

MakeAtmFlux::~MakeAtmFlux()
{

}

/**
 * @brief  		Choose interpolation option and operate
 * @details 		1) BilinearInterpolation
 * @param InterpolationLevel [description]
 */
double MakeAtmFlux::InterpolateFlux(int InterpolationLevel, double Energy, double CosZenith, int AzimuthAngle, int flavor)
{
	if(InterpolationLevel == 1)
	{
		// cout << "Bilinear interpolation" << endl;
		return BilinearInterpolation(Energy, CosZenith, AzimuthAngle, flavor);
	}
	else if(InterpolationLevel == 2)
	{
		cout << "Bicubic interpolation" << endl;
	/*=====================================================================================================================================
	
	Context of bicubicinterpolation 

	 =====================================================================================================================================*/
	}
	else
	{

	}
	return 0;
}

double MakeAtmFlux::BilinearInterpolation(double Energy, double CosZenith, int AzimuthAngle, int flavor)
{
	// We need to consider when the value is in boundary
	// After find the binning of each value, we need to make bilinearinterpolation
	int energyNBin, cosZenithNBin;
	FindProperBin(Energy, CosZenith, &energyNBin, &cosZenithNBin);
	// cout << energyNBin << "\t" << cosZenithNBin << endl;
	// cout << "Energy\t" << Energy << "Cosine\t" << CosZenith << endl;

	// for(int i = 0; i < Azimuth_Length; i++)
	// {
	// 	cout << "azimuth angle : \t" << i * 30 << "\t" <<
	// 	 "energy : \t" << hondaEnergy[cosZenithNBin][i][energyNBin] << "\t" <<
	// 	 "numu flux : \t" << hondaNuMu[cosZenithNBin][i][energyNBin] << "\t" <<
	// 	 "numubar flux : \t" << hondaNuMuBar[cosZenithNBin][i][energyNBin] << "\t" << 
	// 	 "nue flux : \t" << hondaNuE[cosZenithNBin][i][energyNBin] << "\t" <<
	// 	 "nuebar flux : \t" << hondaNuEBar[cosZenithNBin][i][energyNBin] << endl;
	// }

	                    
	/*=====================================================================================================================================
	
	Context of interpolation 

	 =====================================================================================================================================*/
	double eBin1, eBin2; // Energy bin between energy range
	double cBin1, cBin2; // Cosine Zenith angle bin between cosine range
	int aValue = AzimuthAngle; // Azimuth value(example)

	// Maybe we need to set in logarithm?

	if(Energy < hondaEnergy[cosZenithNBin][aValue][0]) Energy = hondaEnergy[cosZenithNBin][aValue][0];
	if(Energy > hondaEnergy[cosZenithNBin][aValue][Contents_Length-1]) Energy = hondaEnergy[cosZenithNBin][aValue][Contents_Length-1];
	
	eBin1 = Energy - hondaEnergy[cosZenithNBin][aValue][energyNBin];
	eBin2 = hondaEnergy[cosZenithNBin][aValue][energyNBin+1] - Energy;

	// cout << "ebin1 : " << eBin1 << endl;
	// cout << "ebin2 : " << eBin2 << endl;


	if(CosZenith > cosineZenith[0]) CosZenith = cosineZenith[0];
	if(CosZenith < cosineZenith[CosZenith_Length - 1]) CosZenith = cosineZenith[CosZenith_Length - 1];

	// cout << CosZenith << endl;
	cBin1 = cosineZenith[cosZenithNBin] - CosZenith;
	cBin2 = CosZenith - cosineZenith[cosZenithNBin + 1];

	// cout << cBin1 << endl;
	// cout << cBin2 << endl;
	
	if(eBin1 < 0 || eBin2 < 0)
	{
		cout << "problem in energy binning!" << endl;
		return -1;
	}
	if(cBin1 < 0 || cBin2 < 0)
	{
		cout << "problem in cosine binning!" << endl;
		return -1;
	}

	double eRatio1, eRatio2;
	double cRatio1, cRatio2;

	eRatio1 = eBin1/(eBin1 + eBin2);
	eRatio2 = eBin2/(eBin1 + eBin2);

	cRatio1 = cBin1/(cBin1 + cBin2);
	cRatio2 = cBin2/(cBin1 + cBin2);
  
  	double fluxNuMu, fluxNuMuBar, fluxNuE, fluxNuEBar;                  

  	// Calculate by bilinear interpolation method
  	fluxNuMu = (eRatio1 * hondaNuMu[cosZenithNBin][aValue][energyNBin + 1] + eRatio2 * hondaNuMu[cosZenithNBin][aValue][energyNBin]) * cRatio2
  				+ (eRatio1 * hondaNuMu[cosZenithNBin+1][aValue][energyNBin + 1] + eRatio2 * hondaNuMu[cosZenithNBin+1][aValue][energyNBin]) * cRatio1;
  	fluxNuMuBar = (eRatio1 * hondaNuMuBar[cosZenithNBin][aValue][energyNBin + 1] + eRatio2 * hondaNuMuBar[cosZenithNBin][aValue][energyNBin]) * cRatio2
  				+ (eRatio1 * hondaNuMuBar[cosZenithNBin+1][aValue][energyNBin + 1] + eRatio2 * hondaNuMuBar[cosZenithNBin+1][aValue][energyNBin]) * cRatio1;
  	fluxNuE = (eRatio1 * hondaNuE[cosZenithNBin][aValue][energyNBin + 1] + eRatio2 * hondaNuE[cosZenithNBin][aValue][energyNBin]) * cRatio2
  				+ (eRatio1 * hondaNuE[cosZenithNBin+1][aValue][energyNBin + 1] + eRatio2 * hondaNuE[cosZenithNBin+1][aValue][energyNBin]) * cRatio1;
  	fluxNuEBar = (eRatio1 * hondaNuEBar[cosZenithNBin][aValue][energyNBin + 1] + eRatio2 * hondaNuEBar[cosZenithNBin][aValue][energyNBin]) * cRatio2
  				+ (eRatio1 * hondaNuEBar[cosZenithNBin+1][aValue][energyNBin + 1] + eRatio2 * hondaNuEBar[cosZenithNBin+1][aValue][energyNBin]) * cRatio1;

  	// cout << "This is certain flux in set energy, cosine \t" <<
  	// 	"NuMu : " << fluxNuMu << "\t"
  	// 	"NuMuBar : " << fluxNuMuBar << "\t"
  	// 	"NuE : " << fluxNuE << "\t"
  	// 	"NuEBar : " << fluxNuEBar << endl;
  	
  	switch(flavor)
  	{
  		case 0:
  			return fluxNuMu;
  		case 1:
	  		return fluxNuMuBar;
  		case 2:
	  		return fluxNuE;
  		case 3:
	  		return fluxNuEBar;
  	}

}


double MakeAtmFlux::GetCrossSection(double energy, int flavor)
{
	switch(flavor)
	{
		case 0:	// NuMu
			return 7.0 * energy * pow(10,-43);
		case 1:  // NuMuBar
			return 3.0 * energy * pow(10,-43);
		case 2:  // NuE
			return 7.0 * energy * pow(10,-43);
		case 3:  // NuEBar
			return 7.0 * energy * pow(10,-43);
	}

}

void MakeAtmFlux::FindProperBin(double energy, double cosZenith, int * energyNBin, int * cosZenithNBin)
{
	int checkEnergyBin, checkCosZenithBin;

	/*=====================================================================================================================================
	
		In here, we need to put how to deal with out of boundary
		To solve this problem, we might need extrapolation

	 =====================================================================================================================================*/
	for(int i = 0; i < Contents_Length - 1; i++) // Find certain range includes certain energy
	{
		if(energy == hondaEnergy[0][0][0])
		{
			checkEnergyBin = 0;
		}
		else if(energy < hondaEnergy[0][0][0])
		{
			checkEnergyBin = 0;
		}
		else if(energy > hondaEnergy[0][0][Contents_Length - 1])
		{
			checkEnergyBin = Contents_Length - 1;
		}
		else if(energy > hondaEnergy[0][0][i] && energy <= hondaEnergy[0][0][i+1])
		{
			checkEnergyBin = i;
		}
	}
	for(int i = 0; i < CosZenith_Length; i++) // Find certain range includes certain energy
	{
		if(cosZenith == cosineZenith[CosZenith_Length - 1])
		{
			checkCosZenithBin = CosZenith_Length - 1;
		}
		else if(cosZenith > cosineZenith[0])
		{
			checkCosZenithBin = 0;
		}
		else if(cosZenith < cosineZenith[CosZenith_Length - 1])
		{
			checkCosZenithBin = CosZenith_Length - 1;
		}
		else if(cosZenith > cosineZenith[i+1] && cosZenith <= cosineZenith[i])
		{
			checkCosZenithBin = i;
		}
	}
	* energyNBin = checkEnergyBin;
	* cosZenithNBin = checkCosZenithBin;
	// cout << "Energy : " << energy << "\t" << checkEnergyBin << endl;
	// cout << "CosZenith : " << cosZenith << "\t" << checkCosZenithBin << endl;
}






#define azNumBin 12
#define czNumBin 20
#define _USE_MATH_DEFINES

int main(void)
{
	/* code */
	string cardFile = "../resources/honda-2015-spl-solmin.d";
	MakeAtmFlux * flux = new MakeAtmFlux(cardFile);
	// flux->Test();
	
	
	
	// flux->InterpolateFlux(1, 0.1, 0.95);
	

	double czStart = 1.0;
	double czEnd = 0.0;
	double czBinWidth;
	czBinWidth = (czStart - czEnd) / double(czNumBin);

	double czBin[czNumBin];
	for(int i = 0; i < czNumBin; i++)
	{
		czBin[i] = czStart - czBinWidth * (double)i;
		cout << czBin[i] << endl;
	}
	flux->InterpolateFlux(1, 0.1, 0.85, 0, 0);

	// unit length to calculate unit area
	double dE = 0.05; // 0.05 GeV
	double dCZ = czBinWidth;
	double dAZ = (2*M_PI)/double(azNumBin);

	cout << dE << "\t" << dCZ << "\t" << dAZ << endl;
	double fluxValue;
	for(int i = 0; i < azNumBin; i++)
	{
		fluxValue = flux->InterpolateFlux(1, 0.15, 0.95, i, 0);
		// cout << fluxValue << endl;
	}
	// EnergyBin[]

	return 0;
}

















