#ifndef LARLITE_SHOWERRECO_CXX
#define LARLITE_SHOWERRECO_CXX

#include "showerReco.h"

namespace larlite {
	bool showerReco::initialize() {
		return true;
	}
	bool showerReco::analyze(storage_manager* storage) {
		auto ev_mcp = storage->get_data<event_mcpart>("largeant");
		if (!ev_mcp){
			std::cerr << "Error: Invalid pointer to MCParticle data" << std::endl;
			exit(1);
		}

		// Loop over all of the particles in the event
		for (auto const& mcp : *ev_mcp){
			// Search for EM shower particles (photons, electrons, positrons)
			if (isShowerPart(mcp)){
				// Add this particle and its ancestors to a shower
				particleChain.clear();
				addToShower(mcp, shower, particleChain, *ev_mcp);
			}
		}	
		return true;
	}
	bool showerReco::finalize() {
		return true;
	}
}

// Returns true if the given MCParticle is a EM shower particle: Photon, Electron or Positron
bool isShowerPart(larlite::mcpart const& _mcp){
	int pdg = _mcp.PdgCode();
	if (pdg == 22 || pdg == 11 || pdg == -11){
		return true;
	}
	else{
		return false;
	}
}

// Returns the index of the vector shower which contains the given MCParticle, and -1 if it has yet to be grouped 
int getShower(larlite::mcpart const& _mcp, std::vector<std::vector<larlite::mcpart>> const& _shower){
	int id = _mcp.TrackId();
	int retVal = -1;
	for (unsigned int i=0; i<_shower.size(); i++){
		for (auto const& p : _shower[i]){
			if (p.TrackId() == id){
				retVal = i;	
			}
		}
	}
	return retVal;
}

// Recursive function that keeps adding particles to the particleChain and advancing to the mother of the current particle
// until it reaches a particle that has been assigned to a shower. At which case it adds the whole particle chain to the shower
void addToShower(larlite::mcpart const& _mcp, std::vector<std::vector<larlite::mcpart>> & _shower, std::vector<larlite::mcpart> & _particleChain, larlite::event_mcpart const& *_ev_mcp){
	// If the particle has not yet been assigned to a shower
	int showerId = getShower(_mcp, _shower);
	if (showerId != -1){
		if (isShowerPart(_mcp)){
			// We are still within the shower
			// Add this particle to the particleChain
			_particleChain.push_back(_mcp);

			// Find the mother of this particle
			int motherId = _mcp.Mother();

			for (auto const& p : _ev_mcp){
				if (p.TrackId == motherId){
					addToShower(p, _shower, _ev_mcp);
					break;
				}
			}
		}
		else{
			// We must reached the end of the shower
			_shower.push_back(_particleChain);
		}
	}
}

#endif
