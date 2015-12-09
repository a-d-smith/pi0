/**
 * \file pi0.h
 *
 * \ingroup pi0
 * 
 * \brief Class def header for a class pi0
 *
 * @author rsjones
 */

/** \addtogroup pi0

    @{*/

#ifndef LARLITE_PI0RECO_H
#define LARLITE_PI0RECO_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH2.h"
#include <vector>
#include <iomanip>
#include <cmath>
#include "TROOT.h"

namespace larlite {
	class pi0Reco : public ana_base{
		public:

		// Default constructor
		pi0Reco(){ _name="pi0Reco"; _fout=0;}
	
		// Default destructor
		virtual ~pi0Reco(){}
	
		virtual bool initialize();
		virtual bool analyze(storage_manager* storage);
		virtual bool finalize();

		std::vector<larlite::mcshower> photons;	
		std::vector<larlite::mcshower> pair1;	
		std::vector<larlite::mcshower> pair2;	
		std::vector<double> allP;
		std::vector<double> allDist;
		std::vector<double> PCutPath;
		std::vector<double> distCutPath;
		double mean, amp, sigma, maxDist;
		TMultiGraph *mg;	
		TH2D   	    *h1;
		TH2D        *h2;
		TGraph      *g1;
		TCanvas     *c1;
	};

	double getClosestDist(TVector3 const& _pis, TVector3 const& _pie, TVector3 const& _pjs, TVector3 const& _pje);
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
