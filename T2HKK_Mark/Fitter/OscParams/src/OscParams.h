#ifndef OscParams_h
#define OscParams_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class OscParams {

  public :

   OscParams(char *cardFile);
   ~OscParams();
   void SetParamValue(int index, double value){ nominal[index] = value; };
   void SetParamError(int index, double err){ error[index] = err; };
   void SetParamType(int index, int typ){ type[index] = typ; };
   void SetParamMinimum(int index, double minimum){ min[index] = minimum; };
   void SetParamMaximum(int index, double maximum){ max[index] = maximum; };
   double GetParamValue(int index){ return nominal[index]; };
   double GetParamError(int index){ return error[index]; };
   int GetParamType(int index){ return type[index]; };
   double GetParamMinimum(int index){ return min[index]; };
   double GetParamMaximum(int index){ return max[index]; };
   int GetNOsc(){ return nOsc; };

   

  private :

   double nominal[6];
   double error[6];
   int type[6];  //0==flat 1==gaussian
   double min[6];
   double max[6];
   int nOsc;

};

#endif

