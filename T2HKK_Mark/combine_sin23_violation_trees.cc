{

  TChain *chain = new TChain("results");

  for(int i=0; i<=100; i++) chain->AddFile(Form("/reno/home/heiwais25/T2HKKSensitivity/combine_source/sin23Violation/mt_bisul/sin23Violation_output_run_%d.root",i));

  TFile *fout = new TFile("sin23_violation_1hk_combine.root","RECREATE");
  TTree *tree = chain->CloneTree();
  tree->Write("results");
  fout->Close();
//
}
  
