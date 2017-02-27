{

  TChain *chain = new TChain("results");
  for(int i=0; i<100; i++) 
	  chain->AddFile(Form("/reno/home/heiwais25/T2HKK_Mark/combine_source/source_%d_run_1638.root",i));

  TFile *fout = new TFile("combine_output/cpViolation_POT_3_sigma_MT_BISUL.root","RECREATE");
  TTree *tree = chain->CloneTree();
  tree->Write("results");
  fout->Close();

}
  
