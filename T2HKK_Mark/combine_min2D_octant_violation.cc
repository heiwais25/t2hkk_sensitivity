{

  TChain *chain = new TChain("results");
  for(int i=0; i<100; i++) 
	  // for(int j=0; j<100; j++)chain->AddFile(Form("/reno/home/heiwais25/T2HKKSensitivity/combine_source/min2D_oct/mt_bisul/min2D_oct_1oct_%d_2oct_%d_run_10000_mt_bisul.root",i,j));
	  for(int j=0; j<100; j++)chain->AddFile(Form("/reno/home/heiwais25/T2HKKSensitivity/combine_source/source_%d_%d_run_2241.root",i,j));

  TFile *fout = new TFile("combine_output/min2D_oct_violation_both_2HK.root","RECREATE");
  TTree *tree = chain->CloneTree();
  tree->Write("results");
  fout->Close();

}
  
