{


  TFile *fil = new TFile("histograms_dcp0_nh_hk_2p5.root");

  char pol[2][10] = {"FHC","RHC"};
  char flavor[6][20] = {"numu","numub","nue","nueb","numuxnue","numubxnueb"};
  char mode[3][20] = {"CCQE","CCnQE","NC"};

  char histnames[5][20] = {"SIG","WS","NC","OTHER","TOTAL"};

  TH1D *event_hists[2][5];

  TH1D *tmp = (TH1D*)fil->Get("pred1Rmu_FHC_numu_CCQE");
  for(int i=0; i<2; i++)
    for(int j=0; j<5; j++){
      event_hists[i][j] = (TH1D*)tmp->Clone(Form("1Rmu_%s_%s",pol[i],histnames[j]));
      event_hists[i][j]->Reset();
    }  

  for(int pi=0; pi<2; pi++)
    for(int fi=0; fi<6; fi++)
      for(int mi=0; mi<3; mi++){
        TH1D *tmp1 = (TH1D*)fil->Get(Form("pred1Rmu_%s_%s_%s",pol[pi],flavor[fi],mode[mi]));
        if( ((pi==0 && fi==0) || (pi==1 && fi==1)) && (mi==0 || mi==1) ) event_hists[pi][0]->Add(tmp1,1.0);
        else if( ((pi==1 && fi==0) || (pi==0 && fi==1)) && (mi==0 || mi==1) ) event_hists[pi][1]->Add(tmp1,1.0);
        else if( mi==2 ) event_hists[pi][2]->Add(tmp1,1.0);
        else event_hists[pi][3]->Add(tmp1,1.0);
        event_hists[pi][4]->Add(tmp1,1.0);
      }   


  //std::cout << event_hists[0][4]->Integral() << std::endl;
  //std::cout << event_hists[1][4]->Integral() << std::endl;

  for(int i=0; i<5; i++) std::cout << histnames[i] << " " <<   event_hists[0][i]->Integral() << std::endl;
  for(int i=0; i<5; i++) std::cout << histnames[i] << " " <<   event_hists[1][i]->Integral() << std::endl;
 
  for(int pi=0; pi<2; pi++) 
    for(int hi=0; hi<5; hi++) 
      for(int bi=1; bi<=event_hists[pi][hi]->GetNbinsX(); bi++) 
        event_hists[pi][hi]->SetBinContent(bi, event_hists[pi][hi]->GetBinContent(bi)*0.1/event_hists[pi][hi]->GetXaxis()->GetBinWidth(bi));

   int colors[5] = {2,4,6,8,7};

   char hist_label[4][50] = {"Signal","Wrong Sign Signal","NC","Other"};

   TCanvas *c1 = new TCanvas("c1","c1",600,500);
   gStyle->SetTitleAlign(13); 
   gPad->SetRightMargin(0.05);
   THStack *hs1 = new THStack("hs1","");
   for(int i=3; i>=0; i--){
      event_hists[0][i]->SetFillColor(colors[i]);
      hs1->Add(event_hists[0][i]);
   }
   hs1->Draw();
   hs1->GetXaxis()->SetRangeUser(0,4);
   hs1->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
   hs1->GetYaxis()->SetTitle("Events/100 MeV");
   hs1->SetTitle("1R#mu, Neutrino Mode, L=1100 km, OAA=2.5#circ");
   hs1->Draw();
   TLegend *leg1 = new TLegend(0.5,0.55,0.85,0.85);
   for(int i=0; i<4; i++) leg1->AddEntry(event_hists[0][i],hist_label[i],"f");
   leg1->SetFillColor(0);
   leg1->SetBorderSize(0.);
   leg1->Draw();
   c1->Print("pred1Rmu_numode_2p5_dcp0_nh.pdf");

   TCanvas *c2 = new TCanvas("c2","c2",600,500);
   gStyle->SetTitleAlign(13);
   gPad->SetRightMargin(0.05);
   THStack *hs2 = new THStack("hs2","");
   for(int i=3; i>=0; i--){
      event_hists[1][i]->SetFillColor(colors[i]);
      hs2->Add(event_hists[1][i]);
   }
   hs2->Draw();
   hs2->GetXaxis()->SetRangeUser(0,4);
   hs2->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
   hs2->GetYaxis()->SetTitle("Events/100 MeV");
   hs2->SetTitle("1R#mu, Antineutrino Mode, L=1100 km, OAA=2.5#circ");
   hs2->Draw();
   TLegend *leg2 = new TLegend(0.5,0.55,0.85,0.85);
   for(int i=0; i<4; i++) leg2->AddEntry(event_hists[1][i],hist_label[i],"f");
   leg2->SetFillColor(0);
   leg2->SetBorderSize(0.);
   leg2->Draw();
   c2->Print("pred1Rmu_anumode_2p5_dcp0_nh.pdf");




}
