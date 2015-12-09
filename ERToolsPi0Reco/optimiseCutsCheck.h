/**
 * \file optimiseCutsCheck.h
 *
 * \ingroup ERToolsPi0Reco
 * 
 * \brief Class def header for a class optimiseCutsCheck
 *
 * @author rsjones
 */

/** \addtogroup ERToolsPi0Reco

    @{*/

#ifndef LARLITE_OPTIMISECUTSCHECK_H
#define LARLITE_OPTIMISECUTSCHECK_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcshower.h"
#include "TLeaf.h"

namespace larlite {
  /**
     \class optimiseCutsCheck
     User custom analysis class made by SHELL_USER_NAME
   */
  class optimiseCutsCheck : public ana_base{
  
  public:

    /// Default constructor
    optimiseCutsCheck(){ _name="optimiseCutsCheck"; _fout=0;}

    /// Default destructor
    virtual ~optimiseCutsCheck(){}

    /** IMPLEMENT in optimiseCutsCheck.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in optimiseCutsCheck.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in optimiseCutsCheck.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();
		//TFile *file;
		/*
		TTree *selectedCut;
		TLeaf *cutLe;
    TLeaf *evNoLe;
    TLeaf *pairNoLe;
    TLeaf *isPi0Le;
		*/

		int evNo, entryNo, cutMin, cutMax, Ncuts;

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
