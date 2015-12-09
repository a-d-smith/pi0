/**
 * \file ERAlgopi0Train.h
 *
 * \ingroup ERToolsPi0Train
 * 
 * \brief Class def header for a class ERAlgopi0Train
 *
 * @author rsjones
 */

/** \addtogroup ERToolsPi0Train

	@{*/

#ifndef ERTOOL_ERALGOPI0TRAIN_H
#define ERTOOL_ERALGOPI0TRAIN_H

#include "ERTool/Base/AlgoBase.h"
#include "TNtuple.h"

namespace ertool {

	/**
	 \class ERAlgopi0Train
	 User custom Algorithm class made by kazuhiro
	 */
	class ERAlgopi0Train : public AlgoBase {
	
		public:
	
		/// Default constructor
		ERAlgopi0Train(const std::string& name="ERAlgopi0Train");
	
		/// Default destructor
		virtual ~ERAlgopi0Train(){};
	
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
			double L, R, angle, E;
			int ID, PDG;
		};

		TNtuple *nt;
		std::vector<TVector3> nhat, r0;
	};
			
}
#endif

/** @} */ // end of doxygen group 
