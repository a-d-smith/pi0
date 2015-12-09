#ifndef LARLITE_PI0TRAIN_H
#define LARLITE_PI0TRAIN_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcshower.h"
#include "TVector3.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include <vector>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <string>

namespace larlite {
  class pi0Train : public ana_base{
  
  public:

    /// Default constructor
    pi0Train(){ _name="pi0Train"; _fout=0;}

    /// Default destructor
    virtual ~pi0Train(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    std::vector<larlite::mcshower> photons;
    TH1D *h1;	
    TCanvas *c1;
    TF1 *f1;
    std::vector<double> allM;

  protected:
    
  };
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
