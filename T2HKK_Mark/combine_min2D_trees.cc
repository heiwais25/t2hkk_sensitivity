{

  TChain *chain = new TChain("results");

  for(int i=0; i<100; i++) 
	  for(int j=0; j<100; j++)chain->AddFile(Form("/reno/home/heiwais25/T2HKKSensitivity/combine_source/min2D/hk/min2D_sin23_0.50_%d_dm23_0.0024_%d_run_gr_10000.root",i,j));

  TFile *fout = new TFile("min2D_sin23_0.50_dm23_0.0024_combine_gr.root","RECREATE");
  TTree *tree = chain->CloneTree();
  tree->Write("results");
  fout->Close();

}
  
