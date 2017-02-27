#define OscParams_cxx
#include "OscParams.h"
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

OscParams::OscParams(char *cardFile){

  nOsc = 6;

  cout << cardFile << std::endl;

  ifstream infile(cardFile);
  string line;
  int iter = 0;
  while (getline(infile, line)){
    if(line.find("#")==0) continue;
  //  std::cout << line << std::endl;
    if(iter>5) {
      std::cout << "Too many lines in oscillation parameter card file" << std::endl;
      continue;
    }
    std::stringstream lineStream(line);
    lineStream >> nominal[iter];
    lineStream >> error[iter];
    lineStream >> type[iter];
    lineStream >> min[iter];
    lineStream >> max[iter];
    std::cout << nominal[iter] << " " << error[iter] << " " << type[iter] << " " << min[iter] << " " << max[iter] << std::endl;
    iter++;
  }

}

OscParams::~OscParams(){

}


