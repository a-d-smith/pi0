#ifndef LARLITE_PI0RECO_CXX
#define LARLITE_PI0RECO_CXX

#include "pi0Reco.h"

namespace larlite {
	bool pi0Reco::initialize() {
		// Parameters of the photon pair invariant mass gaussian
		// These ain't right
		amp = 0.962185;
		sigma = 3.64595;
		mean   = 134.977;

		// x axis: probability of a photon pair coming from a pi0
		// y axis: closest distance of approach
		// the maximum calculated distance is 70.5855
		//mg = new TMultiGraph();
		h1 = new TH2D("h1", "Parameter Space Distribution", 10, 0, 1, 10, 0, 80);
		h2 = new TH2D("h2", "Goodness of Reconstruction Surface", 100, 0, 1, 100, 0, 20);
		//h2 = new TH2D("h2", "Goodness of Reconstruction Surface", 138, 0, 138, 302, 0, 302);
		c1 = new TCanvas("c1", "Parameter Space Distribution", 800, 600);
		
		return true;
	}
  
	bool pi0Reco::analyze(storage_manager* storage) {
		// Load in the shower data
		auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
		
		// Check the loaded data is okay
		if (!ev_mcs){
			std::cerr << "Could not load in the shower data!" << std::endl;
			exit(1);
		}

		photons.clear();

		// Loop over the showes in the event
		for (auto const& mcs : *ev_mcs){
			// If we have a photon shower
			if (mcs.PdgCode() == 22){
				photons.push_back(mcs);	
			}
		}

		int i = 0;
		// Loop of the photons in pairs
		for (auto const& pi : photons){
			TVector3 pis(
					pi.Start().X(),
					pi.Start().Y(),
					pi.Start().Z()
				    );
			TVector3 pie(
					pi.Start().X() + pi.Start().Px(),
					pi.Start().Y() + pi.Start().Py(),
					pi.Start().Z() + pi.Start().Pz()
				    );
			int j = 0;
			for (auto const& pj : photons){
				// Don't compare with self
				if (i > j){
					TVector3 pjs(
						pj.Start().X(),
						pj.Start().Y(),
						pj.Start().Z()
					            );
					TVector3 pje(
						pj.Start().X() + pj.Start().Px(),
						pj.Start().Y() + pj.Start().Py(),
						pj.Start().Z() + pj.Start().Pz()
					            );
					// Calculate distance of closest approach
					double dist = getClosestDist(pis, pie, pjs, pje);
					
					// Calculate the invariant mass of the photon pair
					double Ei = pi.Start().E();
					double Ej = pj.Start().E();
					
					TVector3 piDir = (1/(pie - pis).Mag())*(pie - pis);
					TVector3 pjDir = (1/(pje - pjs).Mag())*(pje - pjs);
					double cosTheta = piDir.Dot(pjDir);
					double M = pow( (2*Ei*Ej*(1-cosTheta)) , 0.5);
					double P = amp*exp(-0.5*pow( (M - mean) /sigma , 2));
					
					h1->Fill(P,dist);

					// Add the photon pair to the list of pairs along with the
					// P-values (probability of originating from a pi0 due to 
					// invariant mass considerations) and distance of closest
					// approach.
					pair1.push_back(pi);					
					pair2.push_back(pj);					
					allP.push_back(P);
					allDist.push_back(dist);

					std::cout << "(i, j) : " << i << ", " << j << std::endl;
					std::cout << "------- M  = " << M << std::endl;
					std::cout << "------- P  = " << P << std::endl;
					std::cout << "---- dist  = " << dist << std::endl;
				}
				j++;			
			}
			i++;
		}

		std::cout << "-----------------------------------------" << std::endl;
		return true;
	}

	bool pi0Reco::finalize() {
		/*
		c1->SetLogz();	
		h1->SetStats(0);
		h1->Draw("colz");
		h1->GetXaxis()->SetTitle("P(originating from #pi^{0})");
		h1->GetYaxis()->SetTitle("Closest distance of approach / cm");
		c1->SaveAs("images/parameterSpaceDistribution.png");
		*/

		// Total number of pairs of photons
		int T = pair1.size();

		int N = 0;
		// Calculate the total number of pi0s
		for (int i=0; i<T; i++){
			if (pair1[i].MotherPdgCode() == 111 && pair2[i].MotherPdgCode() == 111 && pair1[i].MotherTrackID() == pair2[i].MotherTrackID()){
				N++;
			}
		}
		
		bool keepGoing = true;

		// Set the initial guess of the optimum solution
		double PCut = 0.1;
		double distCut = 10;
			
		PCutPath.push_back(PCut);
		distCutPath.push_back(distCut);

		double PCutTry, distCutTry;

		// Define the step size in P and dist
		double Pstep = 0.01;
		double dstep = 0.01;

		// Number of pi0s reconstructed
		int R;
		// Number of non-pi0s reconstructed as pi0s
		int S;

		// Goodness of reconstruction parameter
		double G;
		double GTry;
		int Rbest = 0;
		int Sbest = T-N;
		int iStep, jStep;
		bool firstTry = true;

		while (keepGoing) {
			keepGoing = false;
			iStep = 0;
			jStep = 0;
			// Look around in all directions in the cut-space
			for (int i=-1; i<=1; i+=1){
				for (int j=-1; j<=1; j+=1){
					if (i!=0 && j!=0){
						std::cout << "Looking (" << i << ", " << j << ")" << std::endl;
						PCutTry = PCut + i*Pstep;
						distCutTry = distCut + j*dstep;
				
						// Calculate how many pi0's we reconstructed using the given cut parameters
						R = 0;			
						S = 0;
						for (int i=0; i<T; i++){
							// See if the pair would be reconstructed as a pi0
							if (allP[i] >= PCutTry && allDist[i] <= distCutTry){
								// See if the pair actually came from the same pi0 
								if (pair1[i].MotherPdgCode() == 111 && pair2[i].MotherPdgCode() == 111 && pair1[i].MotherTrackID() == pair2[i].MotherTrackID()){
									R++;
								}
								else{
									S++;
								}
							}
						}
			// Define a "goodness of reconstruction" parameter
			// We have, T pairs of photons
			// 		- of which N are due to pions
			// 			- of which R were reconstructed correctly			-> good
			// 			- of which (N-R) were missed by the reconstruction		-> bad
			// 		- of which (T-N) are not due to pions
			// 			- of which S were incorrectly reconstructed as pions		-> bad
			// 			- of which (T-N-S) were (correctly) not reconstructed as pions 	-> good
			
						double goodThings = R * (T-N-S);
						double normFactor = N * (T-N);
						GTry = goodThings/normFactor;
						std::cout << "R  " << R << std::endl;
						std::cout << "T  " << T << std::endl;
						std::cout << "N  " << N << std::endl;
						std::cout << "S  " << S << std::endl;
						// If all pions were correctly reconstructed    : R = N
						//    and all non-pions were not reconstructed  : S = 0
						//    -> G = 1;
						// If no pions were correctly reconstructed     : R = 0
						//    -> G = 0;
						// If everything was reconstructed as a pion    : R = N, S = T-N
						//    -> G = 0;
	
						if (firstTry){
							std::cout << "First!" << std::endl;
							G = GTry;
							firstTry = false;
						}
						else{
							// Save this step if it improves G
							if (GTry > G){
								G = GTry;
								Rbest = R;
								Sbest = S;
								iStep = i;
								jStep = j;
								keepGoing = true;
								std::cout << "Better!" << std::endl;
							}
						}
						std::cout << "G  " <<  G << std::endl;
						std::cout << "----------" << std::endl;
					}
				}
			}
			if (keepGoing){
				PCut    += iStep * Pstep;
				distCut += jStep * dstep;
				PCutPath.push_back(PCut);
				distCutPath.push_back(distCut);
				std::cout << "Stepped (" << iStep << ", " << jStep << ")" << std::endl;
			}
			else{
				R = Rbest;
				S = Sbest;
			}
		}

		std::cout << "With the cuts: " << std::endl;
		std::cout << "  P >= " << PCut << std::endl;
		std::cout << "  Distance of closest approach <= " << distCut << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl;
		std::cout << R << " / " << N << " pi0s were reconstructed correctly" << std::endl;
		std::cout << S << " / " << T-N << " non-pi0s were incorrectly reconstructed as pi0s" << std::endl;
		std::cout << "Goodness of reconstruction G = " << G << std::endl;

		// Create a surface plot to show the goodness of reco surface
		for (int k=0; k<=100; k++){
			for (int j=0; j<=100; j++){
				// Calculate how many pi0's we reconstructed using the given cut parameters
				R = 0;			
				S = 0;
				double PCut    = k/100;
				double distCut = j/5;
				for (int i=0; i<T; i++){
					// See if the pair would be reconstructed as a pi0
					if (allP[i] >= PCut && allDist[i] <= distCut){
						// See if the pair actually came from the same pi0 
						if (pair1[i].MotherPdgCode() == 111 && pair2[i].MotherPdgCode() == 111 && pair1[i].MotherTrackID() == pair2[i].MotherTrackID()){
							R++;
						}
						else{
							S++;
						}
					}
				}
				//R = k;
				//S = j;
				double goodThings = R * (T-N-S);
				double normFactor = N * (T-N);
				double G = goodThings/normFactor;
				double logG = log(G);
			
				h2->SetBinContent(k,j,logG);
			}		
			std::cout << k << "%" << std::endl;
		}

		double PCutPathArr[PCutPath.size()];
		double distCutPathArr[PCutPath.size()];
		for (int i=0; i<PCutPath.size(); i++){
			PCutPathArr[i] = PCutPath[i];
			distCutPathArr[i] = distCutPath[i];
		}
		g1 = new TGraph(PCutPath.size(), PCutPathArr, distCutPathArr);
		g1->GetXaxis()->SetRangeUser(0,1);
		g1->GetYaxis()->SetRangeUser(0,20);
		h2->SetStats(0);
		h2->Draw("colz");
		g1->Draw("same AL");
		g1->GetXaxis()->SetTitle("P(originating from #pi^{0}) cut");
		g1->GetYaxis()->SetTitle("Closest distance of approach cut / cm");
		c1->SaveAs("images/goodnessOfRecoOverlay.png");
		
		return true;
	}

	double getClosestDist(TVector3 const& _pis, TVector3 const& _pie, TVector3 const& _pjs, TVector3 const& _pje){
		TVector3  u(_pie - _pis);
		TVector3  v(_pje - _pjs);
		TVector3 w0(_pis - _pjs);

		double a = u.Dot(u);
		double b = u.Dot(v);
		double c = v.Dot(v);
		double d = u.Dot(w0);
		double e = v.Dot(w0);
		
		TVector3 distVect(w0 + (1/(a*c - b*b))*( (b*e - c*d)*u - (a*e - b*d)*v ) );
		
		return distVect.Mag();
	}
}
#endif
