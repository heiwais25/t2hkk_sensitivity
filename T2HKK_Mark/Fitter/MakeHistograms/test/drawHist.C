{

  gStyle->SetLineStyleString(5,"[]");
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleSize(0.04);
  gStyle->SetLabelSize(0.04);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetPadBorderMode(0);
//  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  gStyle->SetPalette(1); 



// Raw Probability Plots	
/*
	TH1F *_NuEToNuE3f 	= Bob->Get("NuEToNuE3f");
	TH1F *_NuEToNuMu3f 	= Bob->Get("NuEToNuMu3f");
	TH1F *_NuEToNuTau3f 	= Bob->Get("NuEToNuTau3f");
	TH1F *_NuEToNuX3f 	= Bob->Get("NuEToNuX3f");
	
*/
	TFile * Bob = _file0;
	Bob->cd();
	TFile * Bob1 = _file1;
	Bob1->cd();
	TFile * Bob2 = _file2;
	Bob2->cd();
	TFile * Bob3 = _file3;
	Bob3->cd();
	TFile * Bob4 = _file4;
	Bob4->cd();


	TH1F *_FHC1Re10 	= Bob->Get("FHC1Re");
	TH1F *_RHC1Re10 	= Bob->Get("RHC1Re");
	TH1F *_FHC1Rmu10 	= Bob->Get("FHC1Rmu");
	TH1F *_RHC1Rmu10 	= Bob->Get("RHC1Rmu");

	TH1F *_FHC1Re11 	= Bob1->Get("FHC1Re");
	TH1F *_RHC1Re11 	= Bob1->Get("RHC1Re");
	TH1F *_FHC1Rmu11 	= Bob1->Get("FHC1Rmu");
	TH1F *_RHC1Rmu11 	= Bob1->Get("RHC1Rmu");

	TH1F *_FHC1Re12 	= Bob2->Get("FHC1Re");
	TH1F *_RHC1Re12 	= Bob2->Get("RHC1Re");
	TH1F *_FHC1Rmu12 	= Bob2->Get("FHC1Rmu");
	TH1F *_RHC1Rmu12 	= Bob2->Get("RHC1Rmu");

	TH1F *_FHC1Re09 	= Bob3->Get("FHC1Re");
	TH1F *_RHC1Re09 	= Bob3->Get("RHC1Re");
	TH1F *_FHC1Rmu09 	= Bob3->Get("FHC1Rmu");
	TH1F *_RHC1Rmu09 	= Bob3->Get("RHC1Rmu");

	TH1F *_FHC1Re08 	= Bob4->Get("FHC1Re");
	TH1F *_RHC1Re08 	= Bob4->Get("RHC1Re");
	TH1F *_FHC1Rmu08 	= Bob4->Get("FHC1Rmu");
	TH1F *_RHC1Rmu08 	= Bob4->Get("RHC1Rmu");


/*
	TH1F *_NuTauToNuE3f 	= Bob->Get("NuTauToNuE3f");
	TH1F *_NuTauToNuMu3f 	= Bob->Get("NuTauToNuMu3f");
	TH1F *_NuTauToNuTau3f 	= Bob->Get("NuTauToNuTau3f");
	TH1F *_NuTauToNuX3f 	= Bob->Get("NuTauToNuX3f");
*/
//------------------- End of Raw Probability Plots	


  TCanvas *c4 = new TCanvas("c4","",5,5,800,800);
  c4->SetLeftMargin(0.153226);
  c4->SetRightMargin(0.153226);
  c4->SetBottomMargin(0.153226);
  c4->SetTopMargin(0.103226);
  c4->SetRightMargin(0.053226);
  c4->SetFillColor(kWhite);
  c4->SetLogx(0);

  gStyle->SetPalette(1); 
  gStyle->SetPalette(1); 

	
	_RHC1Re10->GetXaxis()->SetTitleSize(0.03);
	_RHC1Re10->GetXaxis()->SetTitleOffset(2.0);
	_RHC1Re10->GetXaxis()->SetRangeUser(0.0, 3.0);
	_RHC1Re10->GetYaxis()->SetRangeUser(0.0, 40.0);
	// _RHC1Re10->GetYaxis()->SetRangeUser(-1.00, 1.00);
	//_RHC1Re10->GetZaxis()->SetRangeUser(0.00, 0.55);
	_RHC1Re10->GetYaxis()->SetTitleOffset(2.0);
	_RHC1Re10->SetTitle("1Re, Anti Neutrino Mode, L=1100km, OAA=1.5#circ");
	// _RHC1Re10->SetTitleSize(0.04);
	_RHC1Re10->GetXaxis()->SetTitle("Reconstructed Energy(GeV)");
	_RHC1Re10->GetYaxis()->SetTitle("Events/100MeV");
	_RHC1Re10->SetLineColor(1);
	_RHC1Re10->SetLineWidth(2);
	_RHC1Re10->Draw();

	_RHC1Re11->SetLineColor(2);
	_RHC1Re11->SetLineWidth(2);
	// _RHC1Re11->Draw("same");

	_RHC1Re12->SetLineColor(4);
	_RHC1Re12->SetLineWidth(2);
	_RHC1Re12->Draw("same");

	_RHC1Re09->SetLineColor(6);
	_RHC1Re09->SetLineWidth(2);
	// _RHC1Re09->Draw("same");

	_RHC1Re08->SetLineColor(8);
	_RHC1Re08->SetLineWidth(2);
	_RHC1Re08->Draw("same");

	double XLegendStp = 0.65;
	double YLegendStp = 0.72;
	double XLegendWidth = 0.1;
	double YLegendWidth = 0.15;
	TLegend *leg = new TLegend(XLegendStp,YLegendStp,XLegendStp + XLegendWidth    ,YLegendStp + YLegendWidth);

	leg->AddEntry(_RHC1Re08 , "energy scale = 0.8","l");
	leg->AddEntry(_RHC1Re10 , "energy scale = 1.0","l");
	leg->AddEntry(_RHC1Re12 , "energy scale = 1.2","l");
	//leg->AddEntry(gr[1] , "HK + KD(374kt)","l");
	//leg->AddEntry(gr[2] , "HK + KD(561kt)","l");
	//leg->AddEntry(gr[3] , "Step #it{f} #rho / 5year","l");
	//leg->AddEntry(gr[4] , "Const #rho / 10year","l");
	//leg->AddEntry(gr[5] , "Step #it{f} #rho / 10year","l");


	leg->SetTextSize(0.03);
	leg->SetBorderSize(0.);
	leg->SetFillColor(0);
	leg->Draw();


	// _FHC1Re12->Draw("same");


 //  TCanvas *c5 = new TCanvas("c5","",5,5,800,800);
 //  c5->SetLeftMargin(0.153226);
 //  c5->SetRightMargin(0.153226);
 //  c5->SetBottomMargin(0.153226);
 //  c5->SetTopMargin(0.103226);
 //  c5->SetRightMargin(0.053226);
 //  c5->SetFillColor(kWhite);
 //  c5->SetLogx(0);

 //  gStyle->SetPalette(1); 
 //  gStyle->SetPalette(1); 

	
	// _FHC1Re11->GetXaxis()->SetTitleSize(0.03);
	// _FHC1Re11->GetXaxis()->SetTitleOffset(2.0);
	// _FHC1Re11->GetXaxis()->SetRangeUser(0.0, 3.0);
	// // _FHC1Re11->GetYaxis()->SetRangeUser(-1.00, 1.00);
	// //_FHC1Re11->GetZaxis()->SetRangeUser(0.00, 0.55);
	// _FHC1Re11->GetYaxis()->SetTitleOffset(2.0);
	// _FHC1Re11->SetTitle("1Re, Neutrino Mode, L=1100km, OAA=1.5#degree");
	// // _FHC1Re11->SetTitleSize(0.04);
	// _FHC1Re11->GetXaxis()->SetTitle("Reconstructed Energy(GeV)");
	// _FHC1Re11->GetYaxis()->SetTitle("Events/100MeV");
	// _FHC1Re11->Draw();



	// _NuMuToNuMu3f->Draw("same");
	



//   TCanvas *c5 = new TCanvas("c5","",5,5,800,800);
//   c5->SetLeftMargin(0.153226);
//   c5->SetRightMargin(0.153226);
//   c5->SetBottomMargin(0.153226);
//   c5->SetTopMargin(0.153226);
//   c5->SetRightMargin(0.153226);
//   c5->SetFillColor(kWhite);
//   c5->SetLogx(0);
// //  c5->BuildLegend();
// //  c5->SetMarkerStyle(22);
	
// //	leg->AddEntry(Bob,"Histogram","l");

// //	_NuMuToNuMu3f->GetYaxis()->SetRangeUser(-1.00, 1.00);
// 	_NuMuToNuMu3f->GetXaxis()->SetTitleSize(0.03);
// 	_NuMuToNuMu3f->GetXaxis()->SetTitleOffset(2.0);
// 	_NuMuToNuMu3f->GetYaxis()->SetTitleSize(0.03);
// 	_NuMuToNuMu3f->GetYaxis()->SetTitleOffset(2.5);
// //	_NuMuToNuMu3f->SetTitle->("Prob 3 Flavor NuMu->NuMu");
// 	_NuMuToNuMu3f->GetXaxis()->SetTitle("Energy [GeV]");
// 	_NuMuToNuMu3f->GetYaxis()->SetTitle("Probability");
//         _NuMuToNuMu3f->Add(_NuMuToNuE3f, -1.0);   
// 	_NuMuToNuMu3f->SetLineColor(4);
// 	_NuMuToNuMu3f->Draw("cont4z");
// //	_NuMuToNuE3f->SetLineColor(2);

// //	_NuMuToNuE3f->Draw("same");
// //	_NuMuToNuTau3f->Draw("same");
// //        _NuMuToNuTau3f->SetLineColor(6);

// //  	leg->Draw("same");
// //        leg->SetTextAlign(13);




//   TCanvas *c6 = new TCanvas("c6","",5,5,800,800);
//   c6->SetLeftMargin(0.153226);
//   c6->SetRightMargin(0.153226);
//   c6->SetBottomMargin(0.153226);
//   c6->SetTopMargin(0.153226);
//   c6->SetRightMargin(0.153226);
//   c6->SetFillColor(kWhite);
//   c6->SetLogx(1);


// 	_NuMuToNuTau3f->GetYaxis()->SetRangeUser(-1.00, 1.00);
// 	_NuMuToNuTau3f->GetXaxis()->SetTitleSize(0.03);
// 	_NuMuToNuTau3f->GetXaxis()->SetTitleOffset(2.0);
// 	_NuMuToNuTau3f->GetYaxis()->SetTitleSize(0.03);
// 	_NuMuToNuTau3f->GetYaxis()->SetTitleOffset(2.5);
// 	//_NuMuToNuMu3f->SetTitle->("Prob 3 Flavor NuMu->NuTau");
// 	_NuMuToNuTau3f->GetXaxis()->SetTitle("Energy [GeV]");
// 	_NuMuToNuTau3f->GetYaxis()->SetTitle("Cosine Zenith Angle");
// 	_NuMuToNuTau3f->Draw("cont4z");


	
// //------------------- End of Raw Probability Plots	




}


















































