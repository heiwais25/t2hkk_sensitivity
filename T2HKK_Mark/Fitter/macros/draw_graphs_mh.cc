{

  TFile *fi[4];

  //fi[0] = new TFile("2det_hk_total/mh_2det_hk_total_combine.root");
  fi[0] = new TFile("2det_mt_bisul_total/mh_2det_mt_bisul_total_combine.root");
  fi[1] = new TFile("2det_2p5_total/mh_2det_2p5_total_combine.root");
  fi[2] = new TFile("2det_2p0_total/mh_2det_2p0_total_combine.root");
  fi[3] = new TFile("2det_1p5_total/mh_2det_1p5_total_combine.root");

  TGraph *gr[4];

  int colors[4] = {1,2,4,6};

  for(int i=0; i<4; i++){
    TTree *results = (TTree*)fi[i]->Get("results");
    int n = results->Draw("sqrt(DeltaChi2):DeltaCP","Hierarchy==0");
    gr[i] = new TGraph(n,results->GetV2(),results->GetV1());
    gr[i]->SetLineColor(colors[i]);
    gr[i]->SetLineWidth(2);
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  gPad->SetRightMargin(0.05);
  gPad->SetGridy();
  TH1D *tmp = new TH1D("tmp","",100,0,2.0*TMath::Pi());
  tmp->SetStats(false);
  tmp->SetMaximum(12.0);
  tmp->SetMinimum(0.0);
  tmp->SetTitle("True Normal Hierarchy");
  tmp->GetXaxis()->SetTitle("#delta_{cp} (rad.)");
  tmp->GetYaxis()->SetTitle("T_{MH}");
  tmp->Draw();
  for(int i=0; i<4; i++) gr[i]->Draw("L");
  TLegend *leg = new TLegend(0.4,0.65,0.8,0.88);
  //leg->AddEntry(gr[0],"HK#times2","l");
  leg->AddEntry(gr[0],"HK+Mt. Bisul","l");
  leg->AddEntry(gr[1],"HK+KD at 2.5#circ","l");
  leg->AddEntry(gr[2],"HK+KD at 2.0#circ","l");
  leg->AddEntry(gr[3],"HK+KD at 1.5#circ","l");
  leg->SetBorderSize(0.);
  leg->SetFillColor(0);
  leg->Draw();

}
  
