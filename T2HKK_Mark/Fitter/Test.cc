void Test(void)
{
   TH1D *fluxRatios[2][4];
   TH1D *pred1Re[2][6][3];
   TH1D *pred1Rmu[2][6][3];
   TH2D *raw1Re[2][6][3];
   TH2D *raw1Rmu[2][6][3];
   TH2D *raw1ReOsc[2][6][3];
   TH2D *raw1RmuOsc[2][6][3];
   TH2D *pred1Re2D[2][6][3];
   TH2D *pred1Rmu2D[2][6][3];
      //OscProb *oscProb;
   double mcWeights[2][6];
   double massWeight;
   double potWeight[2];
   TFile *fflux;
   TFile *fraw[2][6];
   double baseline;
   double density;
   int detID;

	char * cardfile = "input_hk.card";
	ifstream infile(cardfile);

	string line;
	getline(infile, line);
	std::cout << line.c_str() << std::endl;

	std::cout << "Hello" << std::endl;

	  fflux = new TFile(line.c_str(),"READ");
  //fflux = TDirectoryFile::OpenFile(line.c_str(),"READ");
  //fflux.OpenFile(line.c_str(),"READ");
  std::cout << "test" << std::endl;
  fluxRatios[0][0] = (TH1D*)fflux->Get("weight_numode_numu"); 
  fluxRatios[0][1] = (TH1D*)fflux->Get("weight_numode_numub"); 
  fluxRatios[0][2] = (TH1D*)fflux->Get("weight_numode_nue"); 
  fluxRatios[0][3] = (TH1D*)fflux->Get("weight_numode_nueb"); 
  fluxRatios[1][0] = (TH1D*)fflux->Get("weight_anumode_numu");
  fluxRatios[1][1] = (TH1D*)fflux->Get("weight_anumode_numub");
  fluxRatios[1][2] = (TH1D*)fflux->Get("weight_anumode_nue"); 
  fluxRatios[1][3] = (TH1D*)fflux->Get("weight_anumode_nueb");

  cout << "Load Raw histograms" << std::endl;

  for(int i=0; i<2; i++)
    for(int j=0; j<6; j++){
      getline(infile, line);
      mcWeights[i][j] = atof(line.c_str());
      std::cout << "Weights " << i << " " << j << " " << mcWeights[i][j] << std::endl;
      getline(infile, line);
      fraw[i][j] = new TFile(line.c_str());
      //fraw[i][j].OpenFile(line.c_str(),"READ");
      raw1Re[i][j][0] = (TH2D*)fraw[i][j]->Get("enu_erec_1Re_CCQE");
      raw1Re[i][j][1] = (TH2D*)fraw[i][j]->Get("enu_erec_1Re_CCnQE");
      raw1Re[i][j][2] = (TH2D*)fraw[i][j]->Get("enu_erec_1Re_NC");
      raw1Rmu[i][j][0] = (TH2D*)fraw[i][j]->Get("enu_erec_1Rmu_CCQE");
      raw1Rmu[i][j][1] = (TH2D*)fraw[i][j]->Get("enu_erec_1Rmu_CCnQE");
      raw1Rmu[i][j][2] = (TH2D*)fraw[i][j]->Get("enu_erec_1Rmu_NC");
      raw1ReOsc[i][j][0] = (TH2D*)raw1Re[i][j][0]->Clone(Form("raw1ReOsc_%d%d0",i,j));
      raw1ReOsc[i][j][1] = (TH2D*)raw1Re[i][j][1]->Clone(Form("raw1ReOsc_%d%d1",i,j));
      raw1ReOsc[i][j][2] = (TH2D*)raw1Re[i][j][2]->Clone(Form("raw1ReOsc_%d%d2",i,j));
      raw1RmuOsc[i][j][0] = (TH2D*)raw1Rmu[i][j][0]->Clone(Form("raw1RmuOsc_%d%d0",i,j));
      raw1RmuOsc[i][j][1] = (TH2D*)raw1Rmu[i][j][1]->Clone(Form("raw1RmuOsc_%d%d1",i,j));
      raw1RmuOsc[i][j][2] = (TH2D*)raw1Rmu[i][j][2]->Clone(Form("raw1RmuOsc_%d%d2",i,j));
    }  

  char pol[2][10] = {"FHC","RHC"};
  char flav[6][20] = {"numu","numub","nue","nueb","numuxnue","numubxnueb"}; 
  char imode[3][10] = {"CCQE","CCnQE","NC"};

  for(int i=0; i<2; i++)
    for(int j=0; j<6; j++)
      for(int k=0; k<3; k++){
        TAxis *blah = raw1Re[i][j][k]->GetYaxis();
        pred1Re[i][j][k] = new TH1D(Form("pred1Re_%s_%s_%s",pol[i],flav[j],imode[k]),"",blah->GetNbins(),blah->GetXbins()->GetArray());
        blah = raw1Rmu[i][j][k]->GetYaxis();
        pred1Rmu[i][j][k] = new TH1D(Form("pred1Rmu_%s_%s_%s",pol[i],flav[j],imode[k]),"",blah->GetNbins(),blah->GetXbins()->GetArray());
      }  

  getline(infile, line);
  baseline = atof(line.c_str()); 
  getline(infile, line);
  density = atof(line.c_str()); 

  getline(infile, line);
  massWeight = atof(line.c_str())/22.5;
  getline(infile, line);
  potWeight[0] = 1.04*atof(line.c_str());
  getline(infile, line);
  potWeight[1] = 1.04*atof(line.c_str());
  TAxis *blah = raw1Re[0][0][0]->GetYaxis();
  double Nbins = blah->GetNbins();
  double Xbins = blah->GetXbins();
  // double * XbinsArray = blah->GetXbins()->GetArray();
  // double XbinsArray[blah->GetNbins()];
  std::cout << Nbins << std::endl;
  std::cout << Xbins << std::endl;
  // std::cout << XbinsArray[0] << std::end;
  // 
  
  for(int i=0; i<2; i++)
    for(int j=0; j<6; j++){
      int flavId = j%4;
      for(int k=0; k<3; k++){

        for(int xi=1; xi<=raw1Re[i][j][k]->GetNbinsX()+1; xi++){
          double enu = raw1Re[i][j][k]->GetXaxis()->GetBinCenter(xi);
          int fbin = fluxRatios[i][flavId]->FindBin(enu);
          double fweight = fluxRatios[i][flavId]->GetBinContent(fbin);
          for(int yi=1; yi<=raw1Re[i][j][k]->GetNbinsY()+1; yi++)
            raw1Re[i][j][k]->SetBinContent(xi,yi,raw1Re[i][j][k]->GetBinContent(xi,yi)*fweight*mcWeights[i][j]*massWeight*potWeight[i]);
        }

        for(int xi=1; xi<=raw1Rmu[i][j][k]->GetNbinsX()+1; xi++){
          double enu = raw1Rmu[i][j][k]->GetXaxis()->GetBinCenter(xi);
          int fbin = fluxRatios[i][flavId]->FindBin(enu);
          double fweight = fluxRatios[i][flavId]->GetBinContent(fbin);
          for(int yi=1; yi<=raw1Rmu[i][j][k]->GetNbinsY()+1; yi++)
            raw1Rmu[i][j][k]->SetBinContent(xi,yi,raw1Rmu[i][j][k]->GetBinContent(xi,yi)*fweight*mcWeights[i][j]*massWeight*potWeight[i]);
        }
      }
    }

  for(int xi=1; xi<=raw1Re[0][2][0]->GetNbinsX()+1; xi++){
          double enu = raw1Re[0][2][0]->GetXaxis()->GetBinCenter(xi);
          int fbin = fluxRatios[0][0]->FindBin(enu);
          double fweight = fluxRatios[0][0]->GetBinContent(fbin);
          // double ExpValue = TMath::Log10(enu);
          
          for(int yi=1; yi<=raw1Re[0][2][0]->GetNbinsY()+1; yi++)
          {
            raw1Re[0][2][0]->SetBinContent(xi,yi,raw1Re[0][2][0]->GetBinContent(xi,yi)*fweight*mcWeights[0][2]*massWeight*potWeight[0]);
            double Contents = raw1Re[0][0][0]->GetBinContent(xi, yi);
            std::cout << enu << "\t"  << fbin << "\t" << fweight << "\t" << Contents << endl;
          }
          
        }    
}