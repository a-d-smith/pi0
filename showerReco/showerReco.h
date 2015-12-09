/**
 * \file showerReco.h
 *
 * \ingroup showerReco
 * 
 * \brief Class def header for a class showerReco
 *
 * @author rsjones
 */

/** \addtogroup showerReco

    @{*/

#ifndef LARLITE_SHOWERRECO_H
#define LARLITE_SHOWERRECO_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcpart.h"

namespace larlite {
  /**
     \class showerReco
     User custom analysis class made by SHELL_USER_NAME
   */
  class showerReco : public ana_base{
  
  public:

    /// Default constructor
    showerReco(){ _name="showerReco"; _fout=0;}

    /// Default destructor
    virtual ~showerReco(){}

    /** IMPLEMENT in showerReco.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in showerReco.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in showerReco.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    std::vector<std::vector<larlite::mcpart>> shower;
    std::vector<larlite::mcpart> particleChain;

  protected:
    
  };
}

bool isShowerPart(larlite::mcpart const& _mcp);
int getShower(larlite::mcpart const& _mcp, std::vector<std::vector<larlite::mcpart>> const& _shower);
void addToShower(larlite::mcpart const& _mcp, std::vector<std::vector<larlite::mcpart>> & _shower, std::vector<larlite::mcpart> & _particleChain, larlite::event_mcpart const& *_ev_mcp);

#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
