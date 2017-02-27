#define SystParams_cxx
#include "SystParams.h"
#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "math.h"
#include "stdlib.h"
#include "TH1D.h"

using namespace std;

SystParams::SystParams(char *cardFile){

  if(cardFile==NULL){
    nSysts = 0;
    return;
  }

  cout << cardFile << std::endl;

  ifstream infile(cardFile);
  string line;
  int iter = 0;
  while (getline(infile, line)){
    if(line.find("#")==0) continue;
    iter++;
  }  

  nSysts = iter;

  infile.clear();                 // clear fail and eof bits
  infile.seekg(0, std::ios::beg);

  type = new SYSTTYPE[nSysts];
  value = new double[nSysts];
  error = new double[nSysts];
  nominal = new double[nSysts];
  min = new double[nSysts];
  max = new double[nSysts];
  fluxMode = new int[nSysts];
  xsecMode = new int[nSysts];
  pidType = new int[nSysts];
  polarityMode = new int[nSysts];
  detID = new int[nSysts];

  covMat = new TMatrixDSym(nSysts);
  TVectorD *paramVec = new TVectorD(nSysts);
 
  iter=0; 
  while (getline(infile, line)){
    if(line.find("#")==0) continue;
    std::stringstream lineStream(line);
    int type_tmp;
    lineStream >> value[iter];
    nominal[iter] = value[iter];
    lineStream >> error[iter];
    lineStream >> min[iter];
    lineStream >> max[iter];
    lineStream >> fluxMode[iter];
    lineStream >> xsecMode[iter];
    lineStream >> pidType[iter];
    lineStream >> polarityMode[iter];
    lineStream >> detID[iter];
    lineStream >> type_tmp;
    type[iter] = (SYSTTYPE)type_tmp;
    this->ConvertBinary(fluxMode[iter]);
    this->ConvertBinary(xsecMode[iter]);
    this->ConvertBinary(polarityMode[iter]);
    this->ConvertBinary(detID[iter]);
    for(int i=0; i<nSysts; i++)
      lineStream >> (*covMat)(iter,i);

    (*paramVec)(iter) = nominal[iter];

    std::cout << nominal[iter] << " " << value[iter] << " " << error[iter] << " " << min[iter] << " " << max[iter] << " " << fluxMode[iter] << " " << xsecMode[iter] << " " << pidType[iter] << " " << polarityMode[iter] << " " << detID[iter] << " " << type[iter] << std::endl;
    iter++;
  }   

  for(int i=0; i<nSysts; i++)
    for(int j=0; j<nSysts; j++)
      (*covMat)(i,j) = (*covMat)(i,j)*error[i]*error[j];

  for(int i=0; i<nSysts; i++)
    (*covMat)(i,i) = (*covMat)(i,i)*1.00001;

  throwParams = new ThrowParams(*paramVec,*covMat);

  covMat->Print();
  covMat->Invert();
  covMat->Print();

  paramVec->Delete();
     
}

SystParams::~SystParams(){

  covMat->Delete();
  delete [] nominal; 
  delete [] value;
  delete [] error;
  delete [] max;
  delete [] min;
  delete [] fluxMode;
  delete [] xsecMode;
  delete [] pidType;
  delete [] polarityMode;
  delete [] detID;
  delete [] type;
 

}

double SystParams::GetChi2(){

  double chi2 = 0;
  for(int i=0; i<nSysts; i++)
    for(int j=0; j<nSysts; j++)
      chi2 += (value[i]-nominal[i])*(value[j]-nominal[j])*(*covMat)(i,j);

  return chi2;

}

void SystParams::ConvertBinary(int &decimal){
  //std::cout << std::endl;
  int binary = 0;
  int iter = 1;
  while(decimal>0){
   int bit = decimal%(int)(pow(10,iter)+0.00001)/(int)(pow(10,iter-1)+0.00001);
   binary += (int)(pow(2,iter-1)+0.00001)*bit;
   decimal = decimal-(int)(pow(10,iter-1)+0.00001)*bit;
   iter++;
   //std::cout << bit << std::endl;
   //std::cout << decimal << std::endl;
 }
 decimal = binary;

 return;
}  

void SystParams::ThrowParameters(){   

  std::vector<double> paramValues;
  throwParams->ThrowSet(paramValues);
  for(int i=0; i<nSysts; i++) value[i] = paramValues[i];

  //std::cout << "Throws " << std::endl;
  //for(int i=0; i<nSysts; i++) std::cout << " " << value[i] << std::endl;

  return;

}

void SystParams::ResetParams(){
  
  for(int i=0; i<nSysts; i++) value[i] = nominal[i];

  return;

}
