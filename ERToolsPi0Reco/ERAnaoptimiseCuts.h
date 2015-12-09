/**
 * \file ERAnaoptimiseCuts.h
 *
 * \ingroup ERToolsPi0Reco
 * 
 * \brief Class def header for a class ERAnaoptimiseCuts
 *
 * @author rsjones
 */

/** \addtogroup ERToolsPi0Reco

    @{*/

#ifndef ERTOOL_ERANAOPTIMISECUTS_H
#define ERTOOL_ERANAOPTIMISECUTS_H

#include "ERTool/Base/AnaBase.h"
#include "TNtuple.h"
namespace ertool {

  /**
     \class ERAnaoptimiseCuts
     User custom Analysis class made by kazuhiro
   */
  class ERAnaoptimiseCuts : public AnaBase {
  
  public:

    /// Default constructor
    ERAnaoptimiseCuts(const std::string& name="ERAnaoptimiseCuts");

    /// Default destructor
    virtual ~ERAnaoptimiseCuts(){}

    /// Reset function
    virtual void Reset();

    /// Function to accept fclite::PSet
    void AcceptPSet(const ::fcllite::PSet& cfg);

    /// Called @ before processing the first event sample
    void ProcessBegin();

    /// Function to evaluate input showers and determine a score
    bool Analyze(const EventData &data, const ParticleGraph &ps);

    /// Called after processing the last event sample
    void ProcessEnd(TFile* fout=nullptr);
	
		// Structure to hold data on a shower
		struct showerData{
      TVector3 start, dir, p;
      double E;
      int ID, PDG;
			ertool::NodeID_t nodeID;
    };
		TNtuple *nT;
	
		double cutMin, cutMax;
		int Ncuts, pairNo;

		std::vector<int> Ni, Nm, Nc;
 };
 double PDF(const double &_M, const double &_D);
}
#endif

/** @} */ // end of doxygen group 
