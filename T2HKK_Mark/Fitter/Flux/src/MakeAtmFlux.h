#ifndef MakeAtmFlux_h
#define MakeAtmFlux_h
#define Contents_Length 101
#define CosZenith_Length 20
#define Azimuth_Length 12
#include <string>
using namespace std;
class MakeAtmFlux
{
	public:
		MakeAtmFlux(string fluxFile);
		~MakeAtmFlux();
		double InterpolateFlux(int InterpolationLevel, double Energy, double CosZenith, int AzimuthAngle, int flavor);
		double BilinearInterpolation(double Energy, double CosZenith, int AzimuthAngle, int flavor);
		double GetCrossSection(double energy, int flavor);
		void FindProperBin(double energy, double cosineZenith, int * energyNBin, int * cosZenithNBin);

	private:
		string flavor[4];
		double cosineZenith[20];
		double hondaEnergy[20][12][101];
		double hondaNuMu[20][12][101];
		double hondaNuMuBar[20][12][101];
		double hondaNuE[20][12][101];
		double hondaNuEBar[20][12][101];
};

#endif