#ifndef ERTOOL_ERALGOPI0RECO_CXX
#define ERTOOL_ERALGOPI0RECO_CXX

#include "ERAlgopi0Reco.h"

namespace ertool {

	ERAlgopi0Reco::ERAlgopi0Reco(const std::string& name) : AlgoBase(name){}

	void ERAlgopi0Reco::Reset(){}

	void ERAlgopi0Reco::AcceptPSet(const ::fcllite::PSet& cfg){}

	void ERAlgopi0Reco::ProcessBegin(){
		// Initialise
		Npi0 = 0;
	}

	bool ERAlgopi0Reco::Reconstruct(const EventData &data, ParticleGraph& graph){
		// Runs for each event
	
		// This vector holds all the data we can get on each shower in a nice format
		std::vector<showerData> shower;	

		// Loop over all showers in the event and extract their data
		for (auto const& sh : graph.GetParticleNodes(RecoType_t::kShower)){
			auto const& sh1 = data.Shower(graph.GetParticle(sh).RecoID());
			auto const& part = graph.GetParticle(sh);

			showerData thisShower;

			TVector3 thisStart(sh1.Start()[0], sh1.Start()[1], sh1.Start()[2]);
			TVector3 thisp(part.Momentum()[0], part.Momentum()[1], part.Momentum()[2]);
			TVector3 thisDir = (1/thisp.Mag())*thisp;

			double thisE     = part.Energy();
			int    thisID    = shower.size();
			int    thisPdg   = part.PdgCode(); 

			thisShower.start = thisStart;
			thisShower.dir   = thisDir;
			thisShower.p     = thisp;
			thisShower.E     = thisE;
			thisShower.ID    = thisID;
			thisShower.PDG   = thisPdg;

			shower.push_back(thisShower);
		}

		// Now we have all of the data we need on the showers, here is the plan:
		//
		// Loop over all pairs of photon showers:
		// 		- Calculate their invariant mass and distance of closest approach
		// 		- Apply selection cuts on these values based on the pi0 training file
		// 		- For those selected, calculate the point of closest approach and 
		// 		  choose this as the point where the pi0 decayed.
		//		- Calculate the 4-momentum of pi0 that decayed
		//		- Create a pi0 particle with these values and add it to the particle
		//		  graph as the mother of the photon pair
		//		- Now we can use the pi0 as a particle in later stages of nnbar reco

		// std::cout << "------------------------------" << std::endl;

		
		for (showerData const& s1 : shower){
			for (showerData const& s2 : shower){
				if (s1.ID > s2.ID){
					if (s1.PDG == 22 && s2.PDG == 22){
						// Defining the parameters for the distance of closest approach calculation
						TVector3 w0 = s1.start - s2.start;
						double a = s1.dir.Dot(s1.dir);
						double b = s1.dir.Dot(s2.dir);
						double c = s2.dir.Dot(s2.dir);
						double d = s1.dir.Dot(w0);
						double e = s2.dir.Dot(w0);	

						TVector3 DVect = w0 + ((b*e - c*d)/(a*c - b*b))*s1.dir - ((a*e - b*d)/(a*c - b*b))*s2.dir;
						double D = DVect.Mag();
	
						// Invariant mass calulation
						double M = std::pow(2 * s1.E * s2.E * (1-b), 0.5);
	
						// Apply a cuts on D and M to determine if we should select this pair as 
						// decay products of a pi0
						// From the training file we found (after correcting for detctor edge effects) that the 
						// invariant masses followed a Gaussian distribution with:
						//
						//        mu    = (130.235 +- 0.319) MeV
						//        sigma = ( 19.971 +- 0.394) MeV
						//
						// The true mass of a pi0 is 134.977 MeV for reference. 
						// For now, we will only accept a mass if it is within 1 sigma of mu
						//
						// From the training file, a threshold was decided on D of D < 5 (cm?)

						double Mmu    = 130.235;
						double Msigma = 19.971;
						double Dmax   = 5;

						if (std::abs(M-Mmu) < Msigma && D < Dmax){
							// Reconstruct a pi0
							Npi0++;
						}

					}
				}	
			}	
		}	
		return true;
	}

	void ERAlgopi0Reco::ProcessEnd(TFile* fout){
		// Finalise
		std::cout << "Number of Pi0's reconstructed: " << Npi0 << std::endl;
	}

}

#endif
