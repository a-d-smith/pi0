#ifndef LARLITE_PI0_CXX
#define LARLITE_PI0_CXX

#include "pi0.h"

namespace larlite {
	bool pi0::initialize() {
		dMmax = 900;
		pi0M   = 134.977;
		c1  = new TCanvas("c1","Photon pair invariant mass distribution", 1000, 800);
		h1  = new TH1D("h1","Photon pair invariant mass distribution", 400, 0, 1000);
		sdSpread = new TH1D("sdSpread", "The spread of standard deviations with varying dM", dMmax, 0, dMmax);
		return true;
	}
  
	bool pi0::analyze(storage_manager* storage) {
		// Load in the shower data
		auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
		
		// Check the loaded data is okay
		if (!ev_mcs){
			std::cerr << "Could not load in the shower data!" << std::endl;
			exit(1);
		}

		photons.clear();

		std::cout << "------------------------" << std::endl;
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
					pi.End().X(),
					pi.End().Y(),
					pi.End().Z()
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
						pj.End().X(),
						pj.End().Y(),
						pj.End().Z()
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
					
					h1->Fill(M);
								
					allM.push_back(M);
				}
				j++;			
			}
			i++;
		}

		return true;
	}

	bool pi0::finalize() {
		for (double dM = 1; dM < dMmax; dM += 1){
			// Calculate the mean
			double totM = 0;
			int N = 0;
			for (int i=0; i < allM.size(); i++){
				if (allM[i] <= pi0M+dM && allM[i] >= pi0M-dM){
					int binN = h1->FindBin(allM[i]);
					if (binN != -1){
						if (h1->GetBinContent(binN) > 1){
							totM += allM[i];
							N++;
						}
					}

				}
			}		
			double mean = totM/N;

			// Calcualte the standard deviation
			double X = 0;
			for (int i=0; i < allM.size(); i++){
				if (allM[i] <= pi0M+dM && allM[i] >= pi0M-dM){
					int binN = h1->FindBin(allM[i]);
					if (binN != -1){
						if (h1->GetBinContent(binN) > 1){
							X += pow(allM[i], 2);	
						}
					}
				}
			}
			X /= N;
			double sd = pow(X - pow(mean, 2), 0.5);
			sdSpread->Fill(dM,sd);
			//std::cout << "dM: " << dM << "  Standard deviation: " << sd << std::endl;
		}
		
		sdSpread->GetXaxis()->SetTitle("Cut-off distance from the mean mass, dM, MeV/c^{2}");
		sdSpread->GetYaxis()->SetTitle("The standard deviation of the invariant mass distribution");
		sdSpread->SetStats(false);
		sdSpread->Draw("L");
		//c1->SaveAs("images/standardDevBeforeCut.png");
		/*
		double amplitude = 950;

		std::cout << "Standard deviation: " << sd        << std::endl;
		std::cout << "              Mean: " << mean      << std::endl; 
		std::cout << "         Amplitude: " << amplitude << std::endl;
			
		f1  = new TF1("f1","[0]*exp( -0.5*(pow( (x - [1])/[2] , 2)) )", 0, 280);
		f1->SetParameter(0, amplitude);
		f1->SetParameter(1, mean);
		f1->SetParameter(2, sd);
		h1->SetLineWidth(2);
		h1->SetLineColor(30);
		h1->SetFillColor(9);
		h1->Draw();
		f1->Draw("same");
		//c1->SetLogy();
		c1->SaveAs("images/hopefullyGaussianFitted.png");
		*/	
		return true;
	}


	/*double getDeltaM(vector<double> _Mass, double _p0, double _p1, double _p2, vector<double> _dM){
		double confidence = 0.9;
		double DeltaM = 0;
		int sizeM = 0;
 		if (sizeM <= _Mass.size()){
			double Gaus = [0]*exp( -0.5*(pow( (_Mass[sizeM] - [1])/[2] , 2)) );
			Gaussian.push_back(Gaus);
		}
		sizeM++;
		return DeltaM;
	}*/




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

	/*
	void fitGaussian(TH1D const& _h, std::vector<double> & _param){
		// The param passed to this function must be a vector of three doubles
		// 	_param[0] = amplitude
		// 	_param[1] = mean
		// 	_param[2] = standard deviation
		// The function will then fit a gaussian to the histogram _h and modify
		// the parameters accordingly
		
  		int binTot = _h->GetNbinsX(); // Total number of bins
  		bool keepGoing = true;
  		double dA = 1;     // amplitude step size
  		double dM = 0.001; // mean step size
  		double dS = 0.1;   // SD step size
		double R0, A, M, S, R;

 	 	while (keepGoing){
 	 		R0 = 0;

 			// Calculate R for the current parameters
 			for (int i=0; i<binTot; i++){
 				binC  = _h->GetBinCenter(binTot);//bin center
 				histY = height of the bin
 				fitY  = evaluate current gaussian at binC
 				// Evaluate goodness of fit parameter
 				R0 += pow(histY-fitY, 2);
 			}
 				
 			// Loop in all directions around the current parameters in fit-space
 			keepGoing = false;
 			for (int Ai = -1; Ai <= 1; Ai += 2){
				for (int Mi = -1; Mi <= 1; Mi += 2){
  					for (int Si = -1; Si <= 1; Si += 2){
  						A = _param[0] + Ai*dA; 
  						M = _param[1] + Mi*dM; 
  						S = _param[2] + Si*dS; 
 
 						// Calculate the goodness of fit parameter in the new position in fit-space
 				 		R = 0;
 						for (i=0; i<binTot; i++){
 							binC  = bin center
 							histY = height of the bin
 							fitY  = evaluate current gaussian at binC
 							// Evaluate goodness of fit parameter
 							R += pow(histY-fitY, 2);
 						}
 
 
 						// If a step in this direction would be an improvement then remember it
 						if (R < R0){
 							R0 = R;
 							newA = A;
 							newM = M;
 							newS = S;
 							keepGoing = true;
 						}							
 					}
 			
 				}
 			}

 			// If stepping would be an improvement then do it!
 			if (keepGoing){
				_param[0] = A;
 				_param[1] = M;
 				_param[2] = S;
 			}
 
 	 	}	
 
 */
 
//	void fitGaussian(TH1D const& _h, std::vector<double> & _param){
		// The param passed to this function must be a vector of three doubles
		// 	_param[0] = amplitude
		// 	_param[1] = mean
		// 	_param[2] = standard deviation
		// The function will then fit a gaussian to the histogram _h and modify
		// the parameters accordingly
		
		/*
 * 		binTot = total number of bins
 * 		binMin = minimum bin value
 * 		binMax = maximum bin value
 * 		keepGoing = true;
 * 		dA = amplitude step size
 * 		dM = mean step size
 * 		dS = SD step size
 *	 	while (keepGoing){
 *	 		
 *			// Calculate R for the current parameters
 *	 		R0 = 0;
 *			for (i=0; i<binTot; i++){
 *				binC  = bin center calculated from binMin, binMax and i. Something like binMin + i*(binMax-binMin)/binTot
 *				histY = height of the bin
 *				fitY  = evaluate current gaussian at binC
 *				// Evaluate goodness of fit parameter
 *				R0 += pow(histY-fitY, 2);
 *			}
 *				
 *			// Loop in all directions around the current parameters in fit-space
 * 			keepGoing = false;
 * 			for (Ai = -1; Ai <= 1; Ai += 2){
 * 				for (Mi = -1; Mi <= 1; Mi += 2){
 * 					for (Si = -1; Si <= 1; Si += 2){
 * 						A = _param[0] + Ai*dA; 
 * 						M = _param[1] + Mi*dM; 
 * 						S = _param[2] + Si*dS; 
 *
 *						// Calculate the goodness of fit parameter in the new position in fit-space
 *				 		R = 0;
 *						for (i=0; i<binTot; i++){
 *							binC  = bin center
 *							histY = height of the bin
 *							fitY  = evaluate current gaussian at binC
 *							// Evaluate goodness of fit parameter
 *							R += pow(histY-fitY, 2);
 *						}
 *
 *
 *						// If a step in this direction would be an improvement then remember it
 *						if (R < R0){
 *							R0 = R;
 *							newA = A;
 *							newM = M;
 *							newS = S;
 *							keepGoing = true;
 *						}							
 *					}
 *			
 *				}
 *			}
 *
 * 			// If stepping would be an improvement then do it!
 * 			if (keepGoing){
 *				_param[0] = A;
 *				_param[1] = M;
 *				_param[2] = S;
 *			}
 *
 *	 	}	
 *
 *
 * 		*/
//	}
}
#endif
