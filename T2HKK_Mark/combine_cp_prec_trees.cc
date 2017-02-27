{

  TChain *chain = new TChain("results");

  for(int i=1; i<6; i++) chain->AddFile(Form("/home/jh/Exercise/t2hkk/New/cp_prec_2det_1p5_nuesyst_%d.root",i));

  TFile *fout = new TFile("cp_prec_2det_hk_nosyst_fixdm2_combine.root","RECREATE");
  TTree *tree = chain->CloneTree();
  tree->Write("results");
  fout->Close();

}
  
