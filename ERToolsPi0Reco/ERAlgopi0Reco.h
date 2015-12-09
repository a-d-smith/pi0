/**
 * \file ERAlgopi0Reco.h
 *
 * \ingroup ERToolsPi0Reco
 * 
 * \brief Class def header for a class ERAlgopi0Reco
 *
 * @author rsjones
 */

/** \addtogroup ERToolsPi0Reco

	@{*/

#ifndef ERTOOL_ERALGOPI0RECO_H
#define ERTOOL_ERALGOPI0RECO_H

#include "ERTool/Base/AlgoBase.h"

namespace ertool {

	/**
	 \class ERAlgopi0Reco
	 User custom Algorithm class made by kazuhiro
	 */
	class ERAlgopi0Reco : public AlgoBase {
	
		public:
	
		/// Default constructor
		ERAlgopi0Reco(const std::string& name="ERAlgopi0Reco");
	
		/// Default destructor
		virtual ~ERAlgopi0Reco(){};
	
		/// Reset function
		void Reset();
	
		/// Function to accept fclite::PSet
		void AcceptPSet(const ::fcllite::PSet& cfg);
	
		/// Called @ before processing the first event sample
		void ProcessBegin();
	
		/// Function to evaluate input showers and determine a score
		bool Reconstruct(const EventData &data, ParticleGraph& graph);
	
		/// Called after processing the last event sample
		void ProcessEnd(TFile* fout=nullptr);

		// Structure to hold the cone data on the showers
		struct showerData{
			TVector3 start, dir, p;
			double E;
			int ID, PDG;
		};

		int Npi0;
	};
			
}
#endif

/** @} */ // end of doxygen group 
