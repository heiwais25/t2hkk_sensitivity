{

  TChain *chain = new TChain("results");

  for(int i=1; i<6; i++) chain->AddFile(Form("New/mh_2det_1p5_nuesyst_%d.root",i));

  TFile *fout = new TFile("mh_2det_1p5_nuesyst_combine.root","RECREATE");
  TTree *tree = chain->CloneTree();
  tree->Write("results");
  fout->Close();

}
  
