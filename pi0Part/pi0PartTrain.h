/**
 * \file pi0PartTrain.h
 *
 * \ingroup pi0Part
 * 
 * \brief Class def header for a class pi0PartTrain
 *
 * @author rsjones
 */

/** \addtogroup pi0Part

    @{*/

#ifndef LARLITE_PI0PARTTRAIN_H
#define LARLITE_PI0PARTTRAIN_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mcshower.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TVector3.h"

namespace larlite {
  /**
     \class pi0PartTrain
     User custom analysis class made by SHELL_USER_NAME
   */
  class pi0PartTrain : public ana_base{
  
  public:

    /// Default constructor
    pi0PartTrain(){ _name="pi0PartTrain"; _fout=0;}

    /// Default destructor
    virtual ~pi0PartTrain(){}

    /** IMPLEMENT in pi0PartTrain.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in pi0PartTrain.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in pi0PartTrain.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    std::vector<larlite::mcpart> photons;
    std::vector<double> allM;
    std::vector<double> allMDist;
    std::vector<double> distDotMom;
    std::vector<int> trajIndex;
    TCanvas *c1;
    TH1D *h1;

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
