#ifndef ERTOOL_ERALGOPI0TRAIN_CXX
#define ERTOOL_ERALGOPI0TRAIN_CXX

#include "ERAlgopi0Train.h"

namespace ertool {

	ERAlgopi0Train::ERAlgopi0Train(const std::string& name) : AlgoBase(name){}

	void ERAlgopi0Train::Reset(){}

	void ERAlgopi0Train::AcceptPSet(const ::fcllite::PSet& cfg){}

	void ERAlgopi0Train::ProcessBegin(){
		// Initialise
		
		nt = new TNtuple("pi0Data","pi0Data","M:D:nContained");
						
		// The TPC (for microboone) has 6 bounding planes with normals given by:
		nhat.push_back(TVector3( 0, 0, 1 ));
		nhat.push_back(TVector3( 0, 0,-1 ));
		nhat.push_back(TVector3( 0, 1, 0 ));
		nhat.push_back(TVector3( 0,-1, 0 ));
		nhat.push_back(TVector3( 1, 0, 0 ));
		nhat.push_back(TVector3(-1, 0, 0 ));

		// And points in the plane
		r0.push_back(TVector3( 0 ,   0 , 1036));
		r0.push_back(TVector3( 0 ,   0 ,   0 ));
		r0.push_back(TVector3( 0 ,  115,   0 ));
		r0.push_back(TVector3( 0 , -115,   0 ));
		r0.push_back(TVector3(256,   0 ,   0 ));
		r0.push_back(TVector3( 0 ,   0 ,   0 ));
	}

	bool ERAlgopi0Train::Reconstruct(const EventData &data, ParticleGraph& graph){
		// Runs for each event
	
		// This vector holds all the data we can get on each shower in a nice format
		std::vector<showerData> shower;	

		// Loop over all showers in the event and extract their data
		for (auto const& sh : graph.GetParticleNodes(RecoType_t::kShower)){
			auto const& sh1 = data.Shower(graph.GetParticle(sh).RecoID());
			auto const& part = graph.GetParticle(sh);

			showerData thisShower;

			TVector3 thisStart(sh1.Start()[0], sh1.Start()[1], sh1.Start()[2]);
			TVector3 thisDir(sh1.Dir()[0], sh1.Dir()[1], sh1.Dir()[2]);
			TVector3 thisp(part.Momentum()[0], part.Momentum()[1], part.Momentum()[2]);
			thisDir = (1/thisDir.Mag())*thisDir;

			double thisL     = sh1.Length();
			double thisR     = sh1.Radius();
			double thisAngle = sh1.Angle();
			double thisE     = part.Energy();
			int    thisID    = shower.size();
			int    thisPdg   = part.PdgCode(); 

			thisShower.start = thisStart;
			thisShower.dir   = thisDir;
			thisShower.p     = thisp;
			thisShower.L     = thisL;
			thisShower.R     = thisR;
			thisShower.angle = thisAngle;
			thisShower.E     = thisE;
			thisShower.ID    = thisID;
			thisShower.PDG   = thisPdg;

			shower.push_back(thisShower);
		}

		// Now we have all of the data we can get on the showers, here is the plan:
		//
		// This code aims to train an algorithm which can be used to reconstruct pi0s
		// In order to do this we calculate values of merit for simulated pi0 events
		// We can then apply cuts to these in order to select pi0s from other similar events
		//
		// We can't see a pi0 in a LArTPC but we can see its decay products through the channel:
		// 		pi0 -> gamma + gamma
		//
		// These photons should come from a common vertex (where the pi0 decayed), to test this
		// we calculate the distance of closest approach D for the photons.
		//
		// The invariant mass of the photon pair should ~equal the mass of the pi0
		//
		// 1) Loop over all pairs of showers in the event
		// 2) For each pair we calculate their:
		// 			- Distance of closest approach, D
		// 			- Invariant mass, M
		// 3) Put all this data into a root file so we can choose cut values
		
		/*
		std::cout << "Shower " << s.ID << std::endl;
		std::cout << "	" << "Start : (" << s.start.X() << ", " << s.start.Y() << ", " << s.start.Z() << ")" << std::endl;
		std::cout << "	" << "Dir   : (" << s.dir.X() << ", " << s.dir.Y() << ", " << s.dir.Z() << ")" << std::endl;
		std::cout << "	" << "L     : " << s.L << std::endl;
		std::cout << "	" << "R     : " << s.R << std::endl;
		std::cout << "	" << "angle : " << s.angle << std::endl;
		std::cout << "	" << "dedx  : " << s.dedx << std::endl;
		std::cout << "	" << "E     : " << s.E << std::endl;
		std::cout << "	" << "t     : " << s.t << std::endl;
		*/

		// std::cout << "------------------------------" << std::endl;

		
		for (showerData const& s1 : shower){
			/*
			std::cout << "Shower " << s1.ID << std::endl;
			std::cout << "	" << "Start : (" << s1.start.X() << ", " << s1.start.Y() << ", " << s1.start.Z() << ")" << std::endl;
			std::cout << "	" << "Dir   : (" << s1.dir.X() << ", " << s1.dir.Y() << ", " << s1.dir.Z() << ")" << std::endl;
			std::cout << "	" << "L     : " << s1.L << std::endl;
			std::cout << "	" << "R     : " << s1.R << std::endl;
			std::cout << "	" << "angle : " << s1.angle << std::endl;
			std::cout << "	" << "dedx  : " << s1.dedx << std::endl;
			std::cout << "	" << "E     : " << s1.E << std::endl;
			std::cout << "	" << "t     : " << s1.t << std::endl;
			*/
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
	
						
						// Create a flag for showers which are entirely held within the active volume of the TPC
			   		// nContainted = 0 -> Both showers exceed the boundaries
			   		// nContainted = 1 -> One  shower exceeds the boundaries
			   		// nContainted = 2 -> No   showers exceed the boundaries
						
						double L1, L2;
						bool cont1 = true;
						bool cont2 = true;
	
						for (int i = 0; i < nhat.size(); i++){
							L1 =  ((r0[i] - s1.start).Dot(nhat[i]))/s1.dir.Dot(nhat[i]);
							L2 =  ((r0[i] - s2.start).Dot(nhat[i]))/s2.dir.Dot(nhat[i]);

							if (L1 < s1.L && L1 > 0){ cont1 = false;}
							if (L2 < s2.L && L2 > 0){ cont2 = false;}
						}
						
						if (cont1 && cont2){
							nt->Fill(M, D, 2);
						}
						else if (!cont1 && !cont2){
							nt->Fill(M, D, 0);
						}
						else{
							nt->Fill(M, D, 1);
						}
					}
				}	
			}	
		}	
		return true;
	}

	void ERAlgopi0Train::ProcessEnd(TFile* fout){
		// Finalise

		nt->Write();
	}

}

#endif
