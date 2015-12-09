#ifndef LARLITE_PI0PARTTRAIN_CXX
#define LARLITE_PI0PARTTRAIN_CXX

#include "pi0PartTrain.h"

namespace larlite {
	bool pi0PartTrain::initialize() {
		c1 = new TCanvas("c1", "Photon pair invariant mass distribution: MCTruth", 800, 600);
		h1 = new TH1D("h1","Photon pair invariant mass distribution: MCTruth", 20, 0, 1);
		return true;
	}
  
	bool pi0PartTrain::analyze(storage_manager* storage) {
		// Load in the data on the mcpart
		auto ev_mcp = storage->get_data<event_mcpart>("largeant");

		// Loop over all of the particles in the file
		for (auto const& mcp : *ev_mcp){
			if (mcp.PdgCode() == 22){
				// Particle is a photon
				int motherId = mcp.Mother();
				for (auto const& mcp2 : *ev_mcp){
					if (mcp2.TrackId() == motherId){
						if (mcp2.PdgCode() == 111){
							// Mother was a pi0
							photons.push_back(mcp);
							trajIndex.push_back(mcp.Trajectory().size());	
						}
						break;
					}
				}
			}
		}

		// Loop over the photons in pairs
		int i = 0;
		for (auto const& pi : photons){
			int j = 0;
			for (auto const& pj : photons){
				if (i > j){
					// Get the direction of the ith photon using momenta
					TVector3 piDir(
							pi.Trajectory()[0].Px(),
							pi.Trajectory()[0].Py(),
							pi.Trajectory()[0].Pz()
							);
					TVector3 deltaDistpi(
							pi.Trajectory()[trajIndex.size()].X() - pi.Trajectory()[0].X(),
							pi.Trajectory()[trajIndex.size()].Y() - pi.Trajectory()[0].Y(),
							pi.Trajectory()[trajIndex.size()].Z() - pi.Trajectory()[0].Z()
							);
					piDir = (1/piDir.Mag())*piDir;	
					deltaDistpi = (1/deltaDistpi.Mag())/deltaDistpi;

					// Get the direction of the jth photon using momenta
					TVector3 pjDir(
							pj.Trajectory()[0].Px(),
							pj.Trajectory()[0].Py(),
							pj.Trajectory()[0].Pz()
							);
					TVector3 deltaDistpj(
							pj.Trajectory()[trajIndex.size()].X() - pj.Trajectory()[0].X(),
							pj.Trajectory()[trajIndex.size()].Y() - pj.Trajectory()[0].Y(),
							pj.Trajectory()[trajIndex.size()].Z() - pj.Trajectory()[0].Z()
							);
					pjDir = (1/pjDir.Mag())*pjDir;	
					deltaDistpj = (1/deltaDistpj.Mag())/deltaDistpj;
					
					// Get the energies of the photons	
					double Ei = pi.Trajectory()[0].E();
					double Ej = pj.Trajectory()[0].E();

					// Get the angle between the photons
					double cosTheta = piDir.Dot(pjDir);
					double cosThetaDist = deltaDistpi.Dot(deltaDistpj);

					// Get the invariant mass of the photon pair
					double M = pow( (2*Ei*Ej*(1-cosTheta)), 0.5 );
					double Mdist = pow( (2*Ei*Ej*(1-cosThetaDist)), 0.5 );
					
					allM.push_back(M);
					allMDist.push_back(Mdist);
					distDotMom.push_back(M.Dot(Mdist));
					h1->Fill(distDotMom);	
				}
				j++;
			}
			i++;
		}
		return true;
	}

	bool pi0PartTrain::finalize() {

		double totM;
		for(double const& M : allM){
			totM += M;
		}
		double mean = totM / allM.size();
		
		std::cout << "Mean = " << mean << std::endl;

		h1->Draw();
//		h1->GetXaxis()->SetTitle("Invariant mass - mean / MeV c^{-2}");
//		h1->GetYaxis()->SetTitle("Number of photon pairs");
		c1->SaveAs("images/distDotMom.png");
		return true;
	}

}
#endif
