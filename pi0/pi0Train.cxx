#ifndef LARLITE_PI0TRAIN_CXX
#define LARLITE_PI0TRAIN_CXX

#include "pi0Train.h"

namespace larlite {

  bool pi0Train::initialize() {
    h1  = new TH1D("h1","Photon pair invariant mass distribution: Momentum Zoomed", 100, -5e-12, 5e-12);
    c1 = new TCanvas("c1","Photon pair invariant mass distribution: Momentum Zoomed",800,600);
    f1 = new TF1("f1", "gaus", -5e-12, 5e-12);
    return true;
  }
  
  bool pi0Train::analyze(storage_manager* storage) {
	// Load in shower data
	auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
        if (!ev_mcs){
        	std::cerr << "Could not load in the shower data!" << std::endl;
                exit(1);
        }

	// Find the photons in the event
        photons.clear();
	for (auto const& mcs : *ev_mcs){
 	       if (mcs.PdgCode() == 22){
        	       photons.push_back(mcs);
               }
	}

	// Loop over photons in pairs pairs
	int i = 0;
	for (auto const& pi : photons){
		
                TVector3 piDir(
                              pi.Start().Px(),
                              pi.Start().Py(),
                              pi.Start().Pz()
                            );
                int j = 0;
                for (auto const& pj : photons){
			if (i > j){
				
			        TVector3 pjDir(
	                	              pj.Start().Px(),
			                      pj.Start().Py(),
			                      pj.Start().Pz()
			                    );
				double Ei = pi.Start().E();
	                        double Ej = pj.Start().E();
	
	                        TVector3 piDirHat = (1/piDir.Mag())*piDir;
	                        TVector3 pjDirHat = (1/pjDir.Mag())*pjDir;
	                        double cosTheta = piDirHat.Dot(pjDirHat);
		                double M = pow( (2*Ei*Ej*(1-cosTheta)) , 0.5);

				h1->Fill(M-134.9766);
				allM.push_back(M-134.9766);
			}
			j++;
		}
		i++;
	} 
    return true;
  }

  bool pi0Train::finalize() {
	// Calculate the mean
 	double totM = 0;
	double maxBin = 0;
        int N = 0;
        for (int i=0; i < allM.size(); i++){
              	int binN = h1->FindBin(allM[i]);
                if (binN != -1){
			if (h1->GetBinContent(binN) > maxBin){
				maxBin = h1->GetBinContent(binN);
			}
                        totM += allM[i];
                        N++;
                }
	}
	double mean = totM/N;

	// Calculate the standard deviation
	double X = 0;	
        for (int i=0; i < allM.size(); i++){
		X += pow(allM[i], 2);
	}
	X /= N;
	double sd = pow(X - pow(mean, 2), 0.5);


	// Calculate the amplitude of the gaussian
	double amp = maxBin/N;

	mean += 134.9766;

	std::cout << "              Mean : " << mean << std::endl;
	std::cout << "Standard Deviation : " << sd   << std::endl;
	std::cout << "         Amplitude : " << amp  << std::endl;

	h1->Draw();
	f1->SetParameter(0, amp*N);
	f1->SetParameter(1, mean-134.9766);
	f1->SetParameter(2, sd);
	f1->Draw("same");
	h1->GetXaxis()->SetTitle("(Invariant mass - mean) / MeV/c^{2}");
	h1->GetYaxis()->SetTitle("Number of photon pairs");
	h1->SetStats(kFALSE);
	c1->SaveAs("images/presentations/massDistributionMomentumZoomed.png");
  	return true;
  }

}
#endif
