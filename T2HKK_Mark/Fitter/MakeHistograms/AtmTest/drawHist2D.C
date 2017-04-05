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
	// TFile * Bob1 = _file1;
	// Bob1->cd();
	// TFile * Bob2 = _file2;
	// Bob2->cd();
	// TFile * Bob3 = _file3;
	// Bob3->cd();
	// TFile * Bob4 = _file4;
	// Bob4->cd();

	TH2F * _Event1Re  = Bob->Get("Event2D1Re");
	TH2F * _Event1Rmu = Bob->Get("Event2D1Rmu");

	TCanvas *c4 = new TCanvas("c4","",5,5,800,800);
	c4->SetLeftMargin(0.153226);
	c4->SetRightMargin(0.153226);
	c4->SetBottomMargin(0.153226);
	c4->SetTopMargin(0.153226);
	c4->SetRightMargin(0.153226);
	c4->SetFillColor(kWhite);
	c4->SetLogx(0);

	gStyle->SetPalette(1); 

	
		_Event1Re->GetXaxis()->SetTitleSize(0.03);
		_Event1Re->GetXaxis()->SetTitleOffset(2.0);
		_Event1Re->GetYaxis()->SetTitleSize(0.03);
		_Event1Re->GetYaxis()->SetRangeUser(0.10, 3.00);
		_Event1Re->GetYaxis()->SetTitleOffset(2.0);

		_Event1Re->GetXaxis()->SetTitle("Cosine Zenith Angle");
		_Event1Re->GetYaxis()->SetTitle("Energy [GeV]");
		_Event1Re->SetTitle("#splitline{The number of single electron ring event}{                (without oscillation)}");

		_Event1Re->SetContour(51);
		_Event1Re->Draw("COLZ2");

	TCanvas *c5 = new TCanvas("c5","",5,5,800,800);
	c5->SetLeftMargin(0.153226);
	c5->SetRightMargin(0.153226);
	c5->SetBottomMargin(0.153226);
	c5->SetTopMargin(0.153226);
	c5->SetRightMargin(0.153226);
	c5->SetFillColor(kWhite);
	c5->SetLogx(0);

	// gStyle->SetPalette(1); 
	gStyle->SetPalette(1); 

	
		_Event1Rmu->GetXaxis()->SetTitleSize(0.03);
		_Event1Rmu->GetXaxis()->SetTitleOffset(2.0);
		_Event1Rmu->GetYaxis()->SetTitleSize(0.03);
		_Event1Rmu->GetYaxis()->SetRangeUser(0.10, 3.00);
		_Event1Rmu->GetYaxis()->SetTitleOffset(2.0);

		_Event1Rmu->GetXaxis()->SetTitle("Cosine Zenith Angle");
		_Event1Rmu->GetYaxis()->SetTitle("Energy [GeV]");
		_Event1Rmu->SetTitle("#splitline{The number of single muon ring event}{                (without oscillation)}");
		_Event1Rmu->SetContour(51);

		_Event1Rmu->Draw("COLZ");


	
// //------------------- End of Raw Probability Plots	




}


















































