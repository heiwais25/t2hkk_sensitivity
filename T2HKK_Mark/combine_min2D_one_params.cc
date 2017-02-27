{

  TChain *chain = new TChain("results");

  //for(int i=0; i<100; i++) 
	  for(int j=0; j<300; j++)chain->AddFile(Form("/reno/home/heiwais25/T2HKKSensitivity/combine_source/min2D_new/hk/min2D_new_sin23_0.45_%d_run_300_2hk_reactor.root",j));

  TFile *fout = new TFile("min2D_sin23_0.45_combine_2hk_reactor.root","RECREATE");
  TTree *tree = chain->CloneTree();
  tree->Write("results");
  fout->Close();
}
  
