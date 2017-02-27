#ifndef SystParams_h
#define SystParams_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TMatrixDSym.h"
#include "ThrowParams.h"

enum SYSTTYPE {ERECNORM,ENUNORM,ESCALE,MATTERDENSITY};

class SystParams {

  public :

   SystParams(char *cardFile);
   void ConvertBinary(int &decimal);
   double GetChi2();
   double GetParamValue(int index){ return value[index]; };
   double GetParamNominal(int index){ return nominal[index]; };
   double GetParamError(int index){ return error[index]; };
   int GetFluxMode(int index){ return fluxMode[index]; };
   int GetXsecMode(int index){ return xsecMode[index]; };
   int GetPolarityMode(int index){ return polarityMode[index]; };
   int GetPIDType(int index){ return pidType[index]; };
   int GetDetectorID(int index){ return detID[index]; };
   bool DoDetector(int detid, int index) { return !((detid & detID[index])==0); };
   SYSTTYPE GetParamType(int index){ return type[index]; };
   void SetParamValue(int index, double val){ value[index] = val; };
   double GetMaxEnergy(int index){ return max[index]; };
   double GetMinEnergy(int index){ return min[index]; };
   int GetNSysts(){ return nSysts; };
   void ResetParams();
   void ThrowParameters();
   ~SystParams();

   

  private :

   SYSTTYPE *type;
   double *value;
   double *nominal;
   double *error;
   double *min;
   double *max;
   int *fluxMode;
   int *xsecMode;
   int *polarityMode;
   int *detID;
   int *pidType;
 
   TMatrixDSym *covMat;
   ThrowParams *throwParams;

   int nSysts;
  
};

#endif

