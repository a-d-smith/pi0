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

#ifndef LARLITE_PI0_H
#define LARLITE_PI0_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TH1.h"
#include <vector>
#include <iomanip>
#include <cmath>
#include "TF1.h"
#include "TROOT.h"

namespace larlite {
	class pi0 : public ana_base{
		public:

		// Default constructor
		pi0(){ _name="pi0"; _fout=0;}
	
		// Default destructor
		virtual ~pi0(){}
	
		virtual bool initialize();
		virtual bool analyze(storage_manager* storage);
		virtual bool finalize();

		std::vector<larlite::mcshower> photons;	
		TCanvas *c1;
		TH1D *h1;
		TF1 *f1;
		TH1D *sdSpread;
		//TH1D *sdSpread;
		double distTol, maxM, dMmax, pi0M;
		std::vector<double> allM;
		//std::vector<double> sDev;
		//std::vector<double> deltaM;
		//std::vector<double> Gaussian;
    
	};

	double getClosestDist(TVector3 const& _pis, TVector3 const& _pie, TVector3 const& _pjs, TVector3 const& _pje);
	//double getDeltaM(vector<double>  _Mass, double _p0, double _p1, double _p2, vector<double> _dM);
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
