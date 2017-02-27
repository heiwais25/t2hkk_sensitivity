/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Because there are error in the process of submitting many runs, we need to run certain run number again
This code could take run number from TFile error code
There are several steps to do for fixing easily

1) Copy the error message to "CombineErrorLog"
2) After cheking the length of error message, run this code with root
3) Then run the Fixing Run. In that cshell script, it will operate certain run number

Date 2017-01-10

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/


#define SystParams_cxx
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

void FixingLog(void)
{
    int nLog = 0;
    char * LogFileName = "CombineErrorLog";
    
    //cout << "Test" << std::endl;
    // Checking the file name
    cout << LogFileName << std::endl;

    ifstream infile(LogFileName);
    string line;
    char ErrorLog[500][200];
    int iter = 0;
    while (!infile.eof())
    {
        //if(line.find("#")==0) continue;
        infile.getline(ErrorLog[iter],200);
        //cout << iter << ": th line " << ErrorLog[iter] << endl;
        iter++;
        
        
    }

    int Standard = 91;
    cout << *(ErrorLog[0]+Standard) << endl;
    cout << *(ErrorLog[0]+Standard + 1) << endl;
    cout << *(ErrorLog[0]+Standard + 2) << endl;
    int FirstValue[500], SecondValue[500];
    int FirstSecondGap1 = 3;
    int FirstSecondGap2 = 7;

    ofstream outfile("FixingRun.txt");

    for(int i = 0; i < iter ; i++)
    {
        if(*(ErrorLog[i] + Standard) != '_')
        {
            FirstValue[i] = (*(ErrorLog[i] + Standard - 1) - 48) * 10 + (*(ErrorLog[i] + Standard) - 48);
            if(*(ErrorLog[i] + Standard + FirstSecondGap1) != '_')
            {
                SecondValue[i] = (*(ErrorLog[i] + Standard - 1 + FirstSecondGap1) - 48) * 10 + (*(ErrorLog[i] + Standard + FirstSecondGap1) - 48);
            }
            else
            {
                SecondValue[i] = (*(ErrorLog[i] + Standard - 1 + FirstSecondGap1) - 48);
            }
        }
        else
        {
            FirstValue[i] = *(ErrorLog[i] + Standard - 1) - 48;
            if(*(ErrorLog[i] + Standard - 1 + FirstSecondGap1) != '_')
            {
                SecondValue[i] = (*(ErrorLog[i] + Standard - 2 + FirstSecondGap1) - 48) * 10 + (*(ErrorLog[i] + Standard - 1 + FirstSecondGap1) - 48);
            }
            else
            {
                SecondValue[i] = (*(ErrorLog[i] + Standard - 2 + FirstSecondGap1) - 48);
            }
        }
        cout << "FirstValue : " << FirstValue[i] << " Second : " << SecondValue[i] << endl;
        outfile << FirstValue[i] << " " << SecondValue[i] <<endl;
    }

    outfile.close();

    cout << "The number of line : " << iter << endl;
    nLog = iter;
 

}
