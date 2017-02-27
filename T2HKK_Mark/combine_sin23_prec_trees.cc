{

  TChain *chain = new TChain("results");

  for(int i=0; i<=100; i++) chain->AddFile(Form("/reno/home/heiwais25/T2HKKSensitivity/combine_source/sin23Precision/hk/sin23Precision_output_run_hk_1_0.6_%d.root",i));

  TFile *fout = new TFile("sin23_prec_hk1_0.6_combine.root","RECREATE");
  TTree *tree = chain->CloneTree();
  tree->Write("results");
  fout->Close();

}
  
