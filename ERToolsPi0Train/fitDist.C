int fitDist(){
	// This is root macro which aims to parametrise the M:D distribution for pi0 reconstruction
	//
	// Functional form of parametrisation:
	// 		P(M, D) = A*exp( - (M - mu)^2 / (2*sigma^2) )*exp(-B*(D - delta))
	//

  TFile f("pi0RecoData.root");
  TTree *pi0Data = (TTree*)f.Get("pi0Data");

  TLeaf *MLe = pi0Data->GetLeaf("M");
  TLeaf *DLe = pi0Data->GetLeaf("D");

	TH2D *h1 = new TH2D("h1","M:D", 100, 0, 120, 100, 0, 9*134.977);

  int N = pi0Data->GetEntries();
  for(int i=0; i<N; i++){
      pi0Data->GetEntry(i);
      h1->Fill(DLe->GetValue(0), MLe->GetValue(0));
  }

  TCanvas *c1 = new TCanvas("canv1", "M:D", 1024, 768);
  h1->SetStats(0);
  h1->GetXaxis()->SetTitle("D");
  h1->GetYaxis()->SetTitle("M");

	// M     -> y
	// D     -> x
	// A     -> [0]
	// mu    -> [1]
	// sigma -> [2]
	// B     -> [3]
	// delta -> [4]
	TF2 *f2 = new TF2("f2", "[0]*exp( - ( pow(y-[1], 2) / (2*pow([2], 2)) ) )*exp(-[3]*(x - [4]))", 0, 120, 0, 9*134.977);
	f2->SetParameter(0, 100);
	f2->SetParameter(1, 134.977);
	f2->SetParameter(2, 20);
	f2->SetParameter(3, 1);
	f2->SetParameter(4, 0);
	
	h1->Fit("f2");
	
	h1->Draw("colz");
	f2->Draw("same");

  c1->SaveAs("imgs/jeremyFile/MDdistFitted.png");

	return 0;
}
