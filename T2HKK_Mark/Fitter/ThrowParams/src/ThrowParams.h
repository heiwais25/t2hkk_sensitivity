#ifndef __THROWPARAM_HH__
#define __THROWPARAM_HH__

#include <iostream>
#include <assert.h>
#include <algorithm>
#include <math.h>

#include <TMath.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TVectorT.h>
#include <TDecompChol.h>

#include <vector>
#include <iostream>
#include <TRandom3.h>
#include <TF1.h>

class ThrowParams 
{
 private:
  typedef TVectorT<double> TVectorD;
  typedef TMatrixTSym<double> TMatrixDSym;
  typedef TMatrixT<double> TMatrixD;

  TVectorD    *pvals;
  TMatrixDSym *covar; 
  TMatrixD    *chel_dec;
  int         npars;
  TRandom3    rand; 
  TF1         *gauss;

 public:
  ThrowParams(TVectorD &parms, TMatrixDSym &covm);
  ~ThrowParams();
   /* {
      if(pvals!=NULL)    pvals->Delete();
      if(covar!=NULL)    covar->Delete();
      if(chel_dec!=NULL) chel_dec->Delete();
    };*/
   
  void SetSeed(int seed = 1867) {rand.SetSeed(seed);};
  int GetSize() {return npars;};
  void ThrowSet(std::vector<double> &parms);

 private:
  void CheloskyDecomp(TMatrixD &chel_mat);
  void StdNormRand(double *z);
   
};
#endif
