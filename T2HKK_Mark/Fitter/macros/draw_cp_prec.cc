{

  TFile *fi[4];

  //fi[0] = new TFile("2det_hk_total/cp_prec_2det_hk_total_combine.root");
  fi[3] = new TFile("2det_1p5_total_test/cp_prec_2det_1p5_total_combine.root");
  fi[0] = new TFile("2det_hk_total_test/cp_prec_2det_hk_total_combine.root");
  //fi[0] = new TFile("2det_mt_bisul_total/cp_prec_2det_mt_bisul_total_combine.root");
  fi[2] = new TFile("2det_2p0_total_test/cp_prec_2det_2p0_total_combine.root");
  fi[1] = new TFile("2det_2p5_total/cp_prec_2det_2p5_total_combine.root");

  TGraph *gr_nh[4];
  TGraph *gr_ih[4];

  for(int fiter=0; fiter<4; fiter++){
    gr_nh[fiter] = new TGraph();
    gr_ih[fiter] = new TGraph();
    TTree *results = (TTree*)fi[fiter]->Get("results");
    TGraph *gr;
    results->SetBranchAddress("DeltaChi2Scan",&gr);
    float dcp;
    int hier;
    results->SetBranchAddress("DeltaCP",&dcp);
    results->SetBranchAddress("Hierarchy",&hier);
    
    int ih_iter = 0;
    int nh_iter = 0;
 
    for(int i=0; i<202; i++){
      results->GetEntry(i);
      double prev_chi2 = 500;
      double prev_dcp = -500;
      double first_dcp = -500;
      double sec_dcp = -500;
      for(int j=0; j<=gr->GetN(); j++){
        double dchi2, deltacp;
        gr->GetPoint(j,deltacp,dchi2);
        if(prev_chi2>1.0 && dchi2<=1.0){
          double b = (prev_chi2-prev_dcp/deltacp*dchi2)/(1.0-prev_dcp/deltacp);
          double a = (dchi2-prev_chi2)/(deltacp-prev_dcp);
          //std::cout << deltacp << " " << b << " " << a << std::endl;
          first_dcp = (1.0-b)/a;
          std::cout << prev_dcp << " " <<  deltacp << " " << b << " " << a << " " << first_dcp <<  std::endl;
        } else if(prev_chi2<1.0 && dchi2>=1.0){
          double b = (prev_chi2-prev_dcp/deltacp*dchi2)/(1.0-prev_dcp/deltacp);
          double a = (dchi2-prev_chi2)/(deltacp-prev_dcp);
          sec_dcp = (1.0-b)/a;
          std::cout << prev_dcp << " " <<  deltacp << " " << b << " " << a << " " << sec_dcp <<  std::endl;
        }
        prev_dcp = deltacp;
        prev_chi2 = dchi2;
      }  
      std::cout << fabs(first_dcp-sec_dcp)/2.0 << std::endl << std::endl;
      if(hier ==0){
        gr_nh[fiter]->SetPoint(nh_iter,dcp,fabs(first_dcp-sec_dcp)/2.0*180.0/TMath::Pi());
        nh_iter++;
      } else {
        gr_ih[fiter]->SetPoint(ih_iter,dcp,fabs(first_dcp-sec_dcp)/2.0*180.0/TMath::Pi());
        ih_iter++;
      }
    }
  }  

  std::cout << "Test 1" << std::endl;

  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  gPad->SetRightMargin(0.05);
  TH1D *hist1 = new TH1D("hist1","",100,0,TMath::Pi()*2.0);
  hist1->SetStats(false);
  hist1->SetMinimum(0);
  hist1->SetMaximum(30);
  hist1->GetXaxis()->SetTitle("#delta_{cp} (rad.)");
  hist1->GetYaxis()->SetTitle("#sigma_{#delta_{cp}} (#circ)");
  hist1->SetTitle("True Normal Hierarchy");
  //hist1->GetYaxis()->SetNdivisions(505);
  hist1->Draw();

  gr_nh[0]->SetLineColor(1);
  gr_nh[0]->SetLineWidth(2);
  gr_nh[0]->Draw("C");
  gr_nh[1]->SetLineWidth(2);
  gr_nh[1]->SetLineColor(2);
  //gr_nh[1]->Draw("C");
  gr_nh[2]->SetLineColor(4);
  gr_nh[2]->SetLineWidth(2);
  gr_nh[2]->Draw("C");
  gr_nh[3]->SetLineWidth(2);
  gr_nh[3]->SetLineColor(6);
  gr_nh[3]->Draw("C");

  TLegend *leg = new TLegend(0.3,0.65,0.8,0.85);
  leg->AddEntry(gr_nh[0],"HK#times2","l");
  //leg->AddEntry(gr_nh[0],"HK+Mt. Bisul","l");
  //leg->AddEntry(gr_nh[1],"HK+KD at 2.5#circ","l");
  leg->AddEntry(gr_nh[2],"HK+KD at 2.0#circ","l");
  leg->AddEntry(gr_nh[3],"HK+KD at 1.5#circ","l");
  leg->SetBorderSize(0.);
  leg->SetFillColor(0);
  leg->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",600,500);
  gPad->SetRightMargin(0.05);
  TH1D *hist2 = new TH1D("hist2","",100,0,TMath::Pi()*2.0);
  hist2->SetStats(false);
  hist2->SetMinimum(0);
  hist2->SetMaximum(30);
  hist2->GetXaxis()->SetTitle("#delta_{cp} (rad.)");
  hist2->GetYaxis()->SetTitle("#sigma_{#delta_{cp}} (#circ)");
  hist2->SetTitle("True Inverted Hierarchy");
 // hist2->GetYaxis()->SetNdivisions(505);
  hist2->Draw();

  gr_ih[0]->SetLineColor(1);
  gr_ih[0]->SetLineWidth(2);
  gr_ih[0]->Draw("L");
  gr_ih[1]->SetLineWidth(2);
  gr_ih[1]->SetLineColor(2);
 // gr_ih[1]->Draw("L");
  gr_ih[2]->SetLineColor(4);
  gr_ih[2]->SetLineWidth(2);
  gr_ih[2]->Draw("L");
  gr_ih[3]->SetLineWidth(2);
  gr_ih[3]->SetLineColor(6);
  gr_ih[3]->Draw("L");

  TLegend *leg2 = new TLegend(0.3,0.65,0.8,0.85);
  leg2->AddEntry(gr_ih[0],"HK#times2","l");
  //leg2->AddEntry(gr_ih[0],"HK+Mt. Bisul","l");
  //leg2->AddEntry(gr_ih[1],"HK+KD at 2.5#circ","l");
  leg2->AddEntry(gr_ih[2],"HK+KD at 2.0#circ","l");
  leg2->AddEntry(gr_ih[3],"HK+KD at 1.5#circ","l");
  leg2->SetBorderSize(0.);
  leg2->SetFillColor(0);
  leg2->Draw();


}     

  
