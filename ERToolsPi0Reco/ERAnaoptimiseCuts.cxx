#ifndef ERTOOL_ERANAOPTIMISECUTS_CXX
#define ERTOOL_ERANAOPTIMISECUTS_CXX

#include "ERAnaoptimiseCuts.h"

namespace ertool {

	ERAnaoptimiseCuts::ERAnaoptimiseCuts(const std::string& name) : AnaBase(name)
	{}

	void ERAnaoptimiseCuts::Reset()
	{}

	void ERAnaoptimiseCuts::AcceptPSet(const ::fcllite::PSet& cfg)
	{}

	void ERAnaoptimiseCuts::ProcessBegin(){
		// Setup an NTuple that holds: cut, Nc, Nm, Ni
		nT = new TNtuple("cutData","cutData","cut:Nc:Nm:Ni:pairNo");	
		
		cutMin = 0.0001;
		cutMax = 0.0112;
		Ncuts  = 100;	
		
		Nc.resize(Ncuts,0);							
		Ni.resize(Ncuts,0);
		Nm.resize(Ncuts,0);
		
		pairNo = 0;
	}

	bool ERAnaoptimiseCuts::Analyze(const EventData &data, const ParticleGraph &ps){
		std::cout << "--------------------" << std::endl;
		
		auto const& mc_graph = MCParticleGraph();

		std::vector<showerData> shower;

		// Loop over all showers in the event and add them to the shower vector
		for (auto const& sh : ps.GetParticleNodes(RecoType_t::kShower)){
			auto const& sh1 = data.Shower(ps.GetParticle(sh).RecoID());
			auto const& part = ps.GetParticle(sh);

			showerData thisShower;

			TVector3 thisStart(sh1.Start()[0], sh1.Start()[1], sh1.Start()[2]);
			TVector3 thisp(part.Momentum()[0], part.Momentum()[1], part.Momentum()[2]);
			TVector3 thisDir = (1/thisp.Mag())*thisp;

			double           thisE      = part.Energy();
			int              thisID     = part.ID();
			ertool::NodeID_t thisNodeID = sh;
			int    thisPdg   = part.PdgCode();

			thisShower.start  = thisStart;
			thisShower.dir    = thisDir;
			thisShower.p      = thisp;
			thisShower.E      = thisE;
			thisShower.ID     = thisID;
			thisShower.nodeID = thisNodeID;
			thisShower.PDG    = thisPdg;

			shower.push_back(thisShower);
		}
		// Loop over a range of cut values {
		for (int n=0; n<Ncuts; n++){
			double cut = cutMin + n*(cutMax-cutMin)/Ncuts;
			// Loop over all unique photon shower pairs 
			for (showerData const& s1 : shower){
				for (showerData const& s2 : shower){
					if (s1.ID > s2.ID){
						if (s1.PDG == 22 && s2.PDG == 22){
							// Calcualte their distance of closest approach, D
							TVector3 w0 = s1.start - s2.start;
							double a = s1.dir.Dot(s1.dir);
							double b = s1.dir.Dot(s2.dir);
							double c = s2.dir.Dot(s2.dir);
							double d = s1.dir.Dot(w0);
							double e = s2.dir.Dot(w0);
	
							TVector3 DVect = w0 + ((b*e - c*d)/(a*c - b*b))*s1.dir - ((a*e - b*d)/(a*c - b*b))*s2.dir;
							double D = DVect.Mag();
							
							// Calculate their invariant mass, M
							double M = std::pow(2 * s1.E * s2.E * (1-b), 0.5);
	
							// Calculate the PDF
							double P = PDF(M,D);

							int s1ParentPDG = mc_graph.GetParticle(mc_graph.GetParticle(s1.nodeID).Parent()).PdgCode();
							int s2ParentPDG = mc_graph.GetParticle(mc_graph.GetParticle(s2.nodeID).Parent()).PdgCode();
							bool correct = (s1ParentPDG == 111 && s2ParentPDG == 111);

							// See how well the reconstruction does
							if (P > cut){
								if (correct){
									// Here we have correctly reconstructed a pion
									Nc[n]++;
								}
								else{
									// Here we have incorectly reconstructed a pion
									Ni[n]++;
								}
							}
							else{
								if (correct){
									// Here we have missed a pion
									Nm[n]++;
								}
							}
							pairNo++;
						}
					}
				}
			}
		}
		return true;
	}

	void ERAnaoptimiseCuts::ProcessEnd(TFile* fout){
		// Save the NTuple to fout
		for (int n=0; n<Ncuts; n++){
			double cut = cutMin + n*(cutMax-cutMin)/Ncuts;
			nT->Fill(cut, Nc[n], Nm[n], Ni[n], pairNo);
		}
		nT->Write();

	}

	double PDF(const double &_M, const double &_D){
		// Define the constants for the fitting
		double A     = 1.74356e2;
		double mu    = 1.18181e2;
		double sigma = 3.24905e1;
		double B     = 9.09156e-1;
		double delta = 2.54256e-1;
		
		// The PDF
		double P = A*std::exp( - ( std::pow(_M-mu, 2) / (2*std::pow(sigma, 2)) ) )*std::exp(-B*(_D - delta));
		return P;
	}

}

#endif
