#include "stdlib.h"
#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

double data[4];
//HK nu, HK nu-bar, KD nu, KD nu-bar
double pred[4];
double nom[4];
double factor[4];
TMatrixD *cov;

bool do_det[2];


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    //std::cout << "alpha : " << par[2] << std::endl;
    // Statistical uncertainties part
    double chi2 = 0.;
    for(int i=0; i<4; i++)
    {
        if(!(do_det[0]) && i<2) continue;
        if(!(do_det[1]) && i>=2) continue;

        pred[i] = par[ (i%2==0 ? 1 : 2) ]*nom[i]*(1.0+factor[i]*par[0]);    // Apply uncertainties (par[0] -> fitting in the region of near to 0)
        //pred[i] = par[ i+1 ]*nom[i]*(1.0+factor[i]*par[0]);
        chi2 += pow(data[i]-pred[i],2)/pred[i];
        // std::cout << chi2 << std::endl;
    }

    // 

    // Systematic uncertainties part
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
          chi2 += (par[i+1]-1.0)*(par[j+1]-1.0)*(*cov)(i,j);
        }
    }

    /*for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      chi2 += (par[i+3]-1.0)*(par[j+3]-1.0)*(*cov)(i,j);*/

    f = chi2;

}


void draw_sens()
{

    factor[0] = 1.0;
    factor[1] = -1.0;
    factor[2] = 3.0;
    factor[3] = -3.0;

    cov = new TMatrixD(2,2);

    double sys_fact = pow(5.0,2);

    (*cov)(0,0) = 0.01*0.01*sys_fact;
    (*cov)(1,1) = 0.01*0.01*sys_fact;
    (*cov)(0,1) = 0.0*sys_fact;
    (*cov)(1,0) = 0.0*sys_fact;
    cov->Print();
    cov->Invert();


    double pot_factor = 0.35; 

    nom[0] = 1800.0*pot_factor; // JD Neutrino
    nom[1] = 1800.0*pot_factor; // JD Anti neutrino
    nom[2] = nom[0]/9.0;         // KD Neutrino
    nom[3] = nom[1]/9.0;         // KD Anti neutrino

    double alpha = 0.1;

    TGraph *gr[3];
    for(int i=0; i<3; i++) 
    {
        gr[i] = new TGraph();
    }
    for(int j=0; j<3; j++)
    {
        if(j==0)
        {
            do_det[0] = true;
            do_det[1] = false;
        } 
        else if(j==1)
        {
            do_det[0] = false;
            do_det[1] = true;
        } 
        else 
        {
            do_det[0] = true;
            do_det[1] = true;
        }
        for(int i=0; i<21; i++) // alpha value from 0 to 0.2
        {
            alpha = ((double)i)*0.01;  

            // Applying uncertainties
            data[0] = nom[0]*(1.0+alpha);   
            data[1] = nom[1]*(1.0-alpha);
            data[2] = nom[2]*(1.0+3.0*alpha); 
            data[3] = nom[3]*(1.0-3.0*alpha); 

            TMinuit *gMinuit = new TMinuit(5);  // Number of parameter is 5
            
            gMinuit->Clear(); // Clear the minuit function
            gMinuit->SetFCN(fcn);
            Double_t arglist[10];
            Int_t ierflg = 0;
            // arglist[0] = 1 : Standard
            // arglist[0] = 2 : try to improve minimum (slow)
            arglist[0] = 1;
            gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
            gMinuit->mnparm(0,"alpha",      0.0, 0.001, -1.0,1.0,ierflg);
            gMinuit->mnparm(1,"nunorm",     1.0, 0.001, 0.0,3.0,ierflg);
            gMinuit->mnparm(2,"nubarnorm", 1.0,  0.001, 0.0,3.0,ierflg);
            gMinuit->mnparm(3,"nunorm1",    1.0, 0.001, 0.0,3.0,ierflg);
            gMinuit->mnparm(4,"nubarnorm1", 1.0, 0.001, 0.0,3.0,ierflg);

            gMinuit->FixParameter(0);   // Fix the alpha value
            if(j==0)
            {
                gMinuit->FixParameter(3);
                gMinuit->FixParameter(4);
            }
            else if(j==1)
            {
                gMinuit->FixParameter(1);
                gMinuit->FixParameter(2);   
            }
            else
            {

            }

            //Do the minimization
            arglist[0] = 500000;
            arglist[1] = 0;
            gMinuit->mnexcm("MIGRAD", arglist ,0,ierflg);

            double ln0,ln1,edm,errdef;
            int nvpar,nparx,icstat;
            gMinuit->mnstat(ln0,edm,errdef,nvpar,nparx,icstat);

            std::cout << "Fmin1 : " << ln0 << std::endl;
            std::cout << "FDem1 : " << edm << std::endl;
            std::cout << "ERRDEF1 : " << errdef << std::endl;
            std::cout << "NPARI1 : " << nvpar << std::endl;
            std::cout << "NPARX1 : " << nparx << std::endl;
            std::cout << "ISTAT1 : " << icstat << std::endl;

            gMinuit->Release(0);
            gMinuit->mnexcm("MIGRAD", arglist ,0,ierflg);
            gMinuit->mnstat(ln1,edm,errdef,nvpar,nparx,icstat);

            std::cout << "Fmin2 : " << ln1 << std::endl;
            std::cout << "FDem2 : " << edm << std::endl;
            std::cout << "ERRDEF2 : " << errdef << std::endl;
            std::cout << "NPARI2 : " << nvpar << std::endl;
            std::cout << "NPARX2 : " << nparx << std::endl;
            std::cout << "ISTAT2 : " << icstat << std::endl;

            gr[j]->SetPoint(i,alpha,sqrt(fabs(ln0-ln1+0.000001)));    // a
            gMinuit->Delete();

        }
    } 

    TGraph *gr_sum = new TGraph();
    for(int i=0; i<gr[0]->GetN(); i++)
    {
        double x,y1,y2;
        gr[0]->GetPoint(i,x,y1); 
        gr[1]->GetPoint(i,x,y2); 
        gr_sum->SetPoint(i,x,sqrt(y1*y1+y2*y2));
    }

    TCanvas *c1 = new TCanvas("c1","c1",600,500);
    gPad->SetRightMargin(0.05);
    gPad->SetGridy();
    TH1D *hist = new TH1D("hist","0.35xStastics",100,0.0,0.2);
    hist->SetMinimum(0.0);
    hist->SetMaximum(10.0);
    hist->SetStats(false);
    hist->GetXaxis()->SetTitle("#alpha");
    hist->GetXaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetTitle("Significance (#sigma)");
    hist->Draw();
    gr_sum->SetMarkerStyle(20);
    gr_sum->SetMarkerColor(8);
    gr[0]->SetMarkerStyle(20);
    gr[1]->SetMarkerStyle(20);
    gr[2]->SetMarkerStyle(20);
    gr[0]->SetMarkerColor(2);
    gr[1]->SetMarkerColor(4);
    gr[0]->Draw("P");
    gr[1]->Draw("P");
    gr[2]->Draw("P");
    TLegend *leg = new TLegend(0.12,0.7,0.4,0.85);
    leg->AddEntry(gr[0],"JD Only","p");
    leg->AddEntry(gr[1],"KD Only","p");
    leg->AddEntry(gr[2],"JD+KD","p");
    leg->AddEntry(gr_sum,"JD+KD(Combined)","p");
    leg->SetBorderSize(0.);
    leg->SetFillColor(0);
    leg->Draw();
    gr_sum->Draw("P");

    std::cout << "Japanese Detector: " << std::endl;
    std::cout << "  alpha=0.05 -> " << gr[0]->Eval(0.05) << std::endl;
    std::cout << "  alpha=0.10 -> " << gr[0]->Eval(0.10) << std::endl;
    std::cout << "  alpha=0.20 -> " << gr[0]->Eval(0.20) << std::endl << std::endl;

    std::cout << "Korean Detector: " << std::endl;
    std::cout << "  alpha=0.05 -> " << gr[1]->Eval(0.05) << std::endl;
    std::cout << "  alpha=0.10 -> " << gr[1]->Eval(0.10) << std::endl;
    std::cout << "  alpha=0.20 -> " << gr[1]->Eval(0.20) << std::endl << std::endl;

    std::cout << "Both Detector: " << std::endl;
    std::cout << "  alpha=0.05 -> " << gr[2]->Eval(0.05) << std::endl;
    std::cout << "  alpha=0.10 -> " << gr[2]->Eval(0.10) << std::endl;
    std::cout << "  alpha=0.20 -> " << gr[2]->Eval(0.20) << std::endl << std::endl;


} 


  
