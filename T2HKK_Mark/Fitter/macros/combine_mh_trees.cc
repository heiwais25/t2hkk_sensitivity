{

  TChain *chain = new TChain("results");

  for(int i=0; i<101; i++) chain->AddFile(Form("2det_mt_bisul_total/mh_2det_mt_bisul_total_%d.root",i));

  TFile *fout = new TFile("2det_mt_bisul_total/mh_2det_mt_bisul_total_combine.root","RECREATE");
  TTree *tree = chain->CloneTree();
  tree->Write("results");
  fout->Close();

}
  
